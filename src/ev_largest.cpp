/*
 * ev_ts.cpp
 * Main Program for Eigenvalue computation with SLEPc and PETSc
 * Solves Laplacian in color and spatial space from each timeslice of a
 * configuration in lime
 * via Eigen::Matrix3cd and computes the Eigenpairs using a shell matrix
 * Each Timeslice is HYP-Smeared, the spectrum is mapped to above one and
 * stretched using the Chebyshev-Polynomial T8 of first kind.
 * Eigenvectors of each timeslice are written to one file, Eigenvalues to
 * another
 * Created on: Aug 26, 2013
 * Author: christopher helmes
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
  //__Initialize MPI__
  int mpistat = 0;
  MPI::Init();
  PetscMPIInt world = 0, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //get number of configuration from last argument to main
  int config = atoi( argv[ (argc-1) ] );
  --argc;

  SlepcInitialize(&argc, &argv, (char*)0, NULL);
  // handle input file
  IO* pars = IO::getInstance();
  pars->set_values("check.txt");
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
  int tstart = 0;
  int tend = 0;
  int todo = 0;
  int *tstarts = new int[world];
  int *tends = new int[world];
  int *todos = new int[world];
  if (rank == 0) {
    // calculate chunks minimal chunk size
    int tmp = L0/world;
    // fill array
    for(int i = 0; i < world; ++i)
      todos[i] = tmp;
    // increase the first chunks until all
    // time slices are distributed
    for(int i = 0; i < L0%world; ++i)
      ++todos[i];
    // fill tstarts and tend;
    tstarts[0] = 0;
    tends[0] = todos[0] - 1;
    for(int i = 1; i < world; ++i) {
      tstarts[i] = tends[i-1] + 1;
      tends[i] = tstarts[i] + todos[i] - 1;
    }
  }
  MPI_Scatter(tstarts, 1, MPI_INT, &tstart, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  MPI_Scatter(tends, 1, MPI_INT, &tend, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  MPI_Scatter(todos, 1, MPI_INT, &todo, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  delete [] tstarts;
  delete [] tends;
  delete [] todos;

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
  //hyp-smearing parameters
  const double ALPHA_1 = pars -> get_float("alpha_1");
  const double ALPHA_2 = pars -> get_float("alpha_2");
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
  double* configuration = new double[todo*V_TS];
  //double* configuration = new double[V_4_LIME];
  ierr = read_lime_gauge_field_doubleprec_timeslices(configuration, conf_name,
      L0, L1, L1, L1, tstart, tend);
  ierr = PetscTime(&v2); CHKERRQ(ierr);
  elapsed_time = v2 - v1;
  PetscPrintf(PETSC_COMM_SELF, "%d: init time %f s\n", rank, elapsed_time);
  //loop over timeslices of a configuration
  for(int ts = 0; ts < todo; ++ts) {
    ierr = PetscTime(&v1); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_SELF, "%d: Initializing time slice %d...\n", rank, ts+tstart);
    //--------------------------------------------------------------------------//
    //                              Data input                                  //
    //--------------------------------------------------------------------------//

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
    ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult_Laplacian2D_Largest);
      CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,
        (void(*)())MatMult_Laplacian2D_Largest);
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
    //std::cout << rank << ": successful" << std::endl;  

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
    write_eig_sys_bin("eigenvectors", config, tstart+ts, nev, eigensystem_fix); 
    //recover spectrum of eigenvalues from acceleration
    std::vector<double> evals_save;
    recover_spectrum_large(nconv, evals_accel, evals_save);
    std::cout << rank << ": " << evals_accel.at(0) << " " << evals_save.at(0) <<std::endl;
    write_eigenvalues_bin("eigenvalues", config, tstart+ts, nev, evals_save);
    write_eigenvalues_bin("phases", config, tstart+ts, nev, phase);
  
    //check trace and sum of eigensystem
    vdv = eigensystem.adjoint() * ( eigensystem );
    sum_single = eigensystem.sum();
    sum = vdv.sum(); 
    trc = vdv.trace();
    std::cout << rank << ": tr(V.adj() * V) = " << trc << std::endl;
    std::cout << rank << ": sum(V.adj() * V) = " << sum << std::endl;
    std::cout << rank << ": sum(V) = " << sum_single << std::endl;
    
    //Display user information on files
    printf("%d: %d phase fixed eigenvectors saved successfully \n", rank, nconv);
    printf("%d: %d eigenvalues saved successfully \n", rank, nconv);
    //__Clean up__
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
  } //end loop over timeslices

  ierr = SlepcFinalize();
  MPI::Finalize();
  return 0;
}

