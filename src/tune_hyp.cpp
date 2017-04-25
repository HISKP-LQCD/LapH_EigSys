/*tune_hyp.cpp
 * Program to tune the parameters of the HYP-Smearing used in Smearing the
 * gaugefields
 *
 * Arguments are the ranges of the two paramertes alpha1 and alpha2 as floats
 */
#include <cstdlib>
#include <string>
#include <Eigen/Core>
#include <mpi.h>
#include <slepceps.h>
#include "petsctime.h"
#include "config_utils.h"
#include "eigensystem.h"
#include "io.h"
#include "navigation.h"
#include "par_io.h"
#include "shell_matop.h"
#include "timeslice.h"
#include "recover_spec.h"
#include "read_write.h"

int main(int argc, char **argv) {
  //--------------------------------------------------------------------------//
  //                             Local variables                              //
  //--------------------------------------------------------------------------//
  PetscMPIInt world = 0, rank = 0;

  //get number of configuration from last argument to main
  int config = atoi( argv[ (argc-1) ] );
  --argc;

  SlepcInitialize(&argc, &argv, (char*)0, NULL);
  // handle input file
  IO* pars = IO::getInstance();
  pars->set_values("parameters.txt");
  if(rank == 0) {
    pars->print_summary();
    printf("calculating config %d\n", config);
  }

  char conf_name [200];
  sprintf(conf_name, "%s/conf.%04d", pars->get_path("conf").c_str(), config );

  PetscErrorCode ierr;
  // print info
  PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d, rank = %d\n", world, rank);
  //Time Tracking
  PetscLogDouble v1;
  PetscLogDouble v2;
  PetscLogDouble elapsed_time;
  // start init timer
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_SELF, "%d: Initializing...\n", rank);
  // calculate the timeslices to work on
  // is done on every communicator because I don't
  // know how else to do
  const int L0 = pars -> get_int("LT");

  // local variables
  Mat A;
  EPS eps;
  EPSType type;
  EPSConvergedReason reason;
  //lattice layout from infile
  const int L1 = pars -> get_int("LX");
  PetscInt V3 = pars -> get_int("V3");
  const int V_TS = pars -> get_int("V_TS");
  //calculation parameters from infile
  const PetscInt nev = pars -> get_int("NEV");
  const int MAT_ENTRIES = pars -> get_int("MAT_ENTRIES");
  const double ALPHA_1_i = pars -> get_float("alpha_1_i");
  const double ALPHA_2_i = pars -> get_float("alpha_2_i");
  const double ALPHA_1_f = pars -> get_float("alpha_1_f");
  const double ALPHA_2_f = pars -> get_float("alpha_2_f");
  const int ITER = pars -> get_int("iter");

  //lookup tables
  //Set up navigation
  Nav* lookup = Nav::getInstance();
  lookup -> init();
  Tslice* slice = Tslice::getInstance();
  slice -> Tslice::init();
  //variables holding eigensystem
  PetscInt nconv;
  PetscScalar eigr;
  PetscScalar eigi;
  std::vector<double> evals_accel;
  std::vector<double> phase;
  Vec xr;
  Vec xi;
  //Matrix to hold the entire eigensystem
  Eigen::MatrixXcd eigensystem(MAT_ENTRIES, nev);
  Eigen::MatrixXcd eigensystem_fix(MAT_ENTRIES, nev);
  Eigen::MatrixXcd vdv(nev,nev);
  std::complex<double> sum_single;
  std::complex<double> sum;
  std::complex<double> trc;

  //Initialize memory for configuration
  double* configuration = new double[L0*V_TS];
  const int ts = 10; 
  ierr = read_lime_gauge_field_doubleprec_timeslices(configuration, conf_name,
      L0, L1, L1, L1, 0, L0 );
  ierr = PetscTime(&v2); CHKERRQ(ierr);
  elapsed_time = v2 - v1;
  PetscPrintf(PETSC_COMM_SELF, "%d: init time %f s\n", rank, elapsed_time);
  //loop over timeslices of a configuration
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  //--------------------------------------------------------------------------//
  //                              Data input                                  //
  //--------------------------------------------------------------------------//
  
  std::vector<std::array<double, 3> > lowest_ev;
  //loop over hyp-smearing parameters, leaving number of iterations constant
  for (double a1 = ALPHA_1_i; a1 <= ALPHA_1_f; a1 += 0.05){
    for (double a2 = ALPHA_2_i; a2 <= ALPHA_2_f; a2 += 0.05){
      std::array<double, 3> meas;
      meas.at(0) = a1;
      meas.at(1) = a2;

      const double ALPHA_1 = a1;
      const double ALPHA_2 = a2;
      //Time Slice of Configuration
      double* timeslice = configuration + (ts*V_TS);

      //Write Timeslice in Eigen Array                                                  
      slice -> map_timeslice_to_eigen(timeslice);
      //Gauge_matrices
      //Gaugetrafo of timeslice
      //slice -> transform_ts(gauge);
      //save transformed timeslice
      //write_link_matrices_ts("ts_gauged_000.1300");
      //Apply Smearing algorithm to timeslice ts
      slice -> smearing_hyp(ALPHA_1, ALPHA_2, ITER);
      //__Define Action of Laplacian in Color and spatial space on vector
      //std::cout << rank << ": Try to create Shell Matrix..." << std::endl;
      ierr = MatCreateShell(PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,MAT_ENTRIES,MAT_ENTRIES,&V3,&A);
        CHKERRQ(ierr);
      //std::cout << rank << ": done" << std::endl;
      //std::cout << rank << ": Try to set operations..." << std::endl;
      ierr = MatSetFromOptions(A);
        CHKERRQ(ierr);
      ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult_Laplacian2D);
        CHKERRQ(ierr);
      ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,
          (void(*)())MatMult_Laplacian2D);
        CHKERRQ(ierr);
      //std::cout << rank << ": accomplished" << std::endl;
      ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,
          (void(*)())MatGetDiagonal_Laplacian2D);
        CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_SELF, "%d: set operations.\n", rank);
      //std::cout << rank << ": Matrix creation..." << std::endl;
        CHKERRQ(ierr);
      ierr = MatSetOption(A, MAT_HERMITIAN, PETSC_TRUE);
        CHKERRQ(ierr);
      ierr = MatSetUp(A);
        CHKERRQ(ierr);
      ierr = MatCreateVecs(A,NULL,&xr);
        CHKERRQ(ierr);
      ierr = MatCreateVecs(A,NULL,&xi);
        CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_SELF, "%d: matrix creation.\n", rank);

      //--------------------------------------------------------------------------//
      //                 Context creation & Options setting                       //
      //--------------------------------------------------------------------------//

      //Create eigensolver context
      ierr = EPSCreate(PETSC_COMM_SELF, &eps);
        CHKERRQ(ierr);

      //Associate A with eps
      ierr = EPSSetOperators(eps, A, NULL);
        CHKERRQ(ierr);
      ierr = EPSSetProblemType(eps, EPS_HEP);
        CHKERRQ(ierr);

      //Set solver parameters at runtime
      //default nev: 200
      ierr = EPSSetDimensions(eps, nev, PETSC_DECIDE, PETSC_DECIDE);
        CHKERRQ(ierr);
      //
      ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE);
        CHKERRQ(ierr);
      //Ask for command line options
      ierr = EPSSetFromOptions(eps);
        CHKERRQ(ierr);

      ierr = PetscTime(&v2); CHKERRQ(ierr);
      elapsed_time = v2 - v1;
      PetscPrintf(PETSC_COMM_SELF, "%d: timing %f s\n", rank, elapsed_time);
      //std::cout << rank << ": solving time " << elapsed_time << std::endl;

      //--------------------------------------------------------------------------//
      //                         Solve Ax = kx                                    //
      //--------------------------------------------------------------------------//

      ierr = PetscTime(&v1); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_SELF, "%d: Start solving T(B) * x = y\n", rank);
      //std::cout << rank << ": Start solving T(B) * x = y" << std::endl;
      ierr = EPSSolve(eps); CHKERRQ(ierr);
      ierr = PetscTime(&v2); CHKERRQ(ierr);
      elapsed_time = v2 - v1;
      PetscPrintf(PETSC_COMM_SELF, "%d: solving time %f s\n", rank, elapsed_time);
      //std::cout << rank << ": solving time " << elapsed_time << std::endl;

      //--------------------------------------------------------------------------//
      //                         Handle solution                                  //
      //--------------------------------------------------------------------------//  
      ierr = EPSGetConvergedReason(eps,&reason);
      CHKERRQ(ierr);
      ierr = EPSGetType(eps,&type);
      CHKERRQ(ierr);
      ierr = EPSGetConverged(eps, &nconv);
      CHKERRQ(ierr);
      //__Retrieve solution__
      ierr = PetscPrintf(PETSC_COMM_SELF,"%D: Number of converged eigenvalues: %D\n", rank, nconv);
      CHKERRQ(ierr);
      ierr=PetscPrintf(PETSC_COMM_SELF, "%D: Convergence reason of solver: %D\n", rank, reason);
      CHKERRQ(ierr);

      nconv = nev;
      //Write nev Eigenpairs to evectors and evalues
      evals_accel.resize(nconv);
      for ( PetscInt i = 0; i < nconv; ++i ) {
        ierr = EPSGetEigenpair(eps,i,&eigr,&eigi,xr,xi);
        CHKERRQ(ierr); 
        PetscScalar* ptr_evec;
        //double* ptr_eval;
        double eval = PetscRealPart(eigr);
        //ptr_eval = &eval;
        ierr = VecGetArray(xr, &ptr_evec);
        CHKERRQ(ierr);
      
        //write each eigenvector to eigensystem for check
        eigensystem.col(i) = Eigen::Map<Eigen::VectorXcd, 0 >(ptr_evec, MAT_ENTRIES);
        //fwrite(ptr_eval ,  sizeof(double) , 1, evalues );
        evals_accel.at(i) = eval;
        ierr = VecRestoreArray(xr, &ptr_evec);
        CHKERRQ(ierr);
      }
      fix_phase(eigensystem, eigensystem_fix, phase);
      //recover spectrum of eigenvalues from acceleration
      std::vector<double> evals_save;
      recover_spectrum(nconv, evals_accel, evals_save);
      std::cout << rank << ": " << evals_accel.at(0) << " " << evals_save.at(0) <<std::endl;    
      // identify eigenvalues and phases by alpha 1 and alpha 2 value
      char ev_prefix [200];
      char phs_prefix [200];
      sprintf(ev_prefix,"eigenvalues_hyp-%d_a1-%.2f_a2-%.2f",ITER,a1,a2);
      sprintf(phs_prefix,"phases_hyp-%d_a1-%.2f_a2-%.2f",ITER,a1,a2);
      write_eigenvalues_bin(ev_prefix, config, ts, nev, evals_save);
      write_eigenvalues_bin(phs_prefix, config, ts, nev, phase);
      meas.at(2) = evals_save.at(0);
      lowest_ev.push_back(meas);
      //check trace and sum of eigensystem
      vdv = eigensystem.adjoint() * ( eigensystem );
      sum_single = eigensystem.sum();
      sum = vdv.sum(); 
      trc = vdv.trace();
      std::cout << rank << ": tr(V.adj() * V) = " << trc << std::endl;
      std::cout << rank << ": sum(V.adj() * V) = " << sum << std::endl;
      std::cout << rank << ": sum(V) = " << sum_single << std::endl;
      
      //Display user information on files
      //printf("%d: %d phase fixed eigenvectors saved successfully \n", rank, nconv);
      printf("%d: %d eigenvalues saved successfully \n", rank, nconv);
    }
  }//end loops over alphas
  // print out lowest eigenvalues and alphas
  std::cout << "#alpha1\talpha2\tlowest_ev" << std::endl;
  for(std::array<double, 3> it : lowest_ev){
    std::cout << it.at(0) << "\t" << it.at(1) << "\t" << it.at(2) << std::endl; 
  }
  //__Clean up__
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}

