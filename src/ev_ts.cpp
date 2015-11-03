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
//__Global Declarations__
//lookup tables
//int up_3d[V3][3],down_3d[V3][3];
//Eigen Array
//Eigen::Matrix3cd **eigen_timeslice = new Eigen::Matrix3cd *[V3];

int main(int argc, char **argv) {
  //--------------------------------------------------------------------------//
  //                             Local variables                              //
  //--------------------------------------------------------------------------//
  //Eigen::initParallel();
  //Eigen::setNbThreads(6);
  //__Initialize MPI__
  int mpistat = 0;
  MPI::Init();

  std::cout << "Values from parameters.txt set" << std::endl; 
  Mat A;              
  EPS eps;             
  EPSType type;
  EPSConvergedReason reason;
  //Handling infile
  IO* pars = IO::getInstance();
  pars -> set_values("parameters.txt");
  pars -> print_summary();
  //Set up navigation
  Nav* lookup = Nav::getInstance();
  lookup -> init();
  //in and outpaths
  std::string GAUGE_FIELDS = pars -> get_path("conf");
  //lattice layout from infile
  const int L0 = pars -> get_int("LT");
  const int L1 = pars -> get_int("LX");
  const int L2 = pars -> get_int("LY");
  const int L3 = pars -> get_int("LZ");
  const int V3 = pars -> get_int("V3");
  const int V_TS = pars -> get_int("V_TS");
  //calculation parameters from infile
  const int NEV = pars -> get_int("NEV");
  const int V_4_LIME = pars -> get_int("V4_LIME");
  const int MAT_ENTRIES = pars -> get_int("MAT_ENTRIES");
  //chebyshev parameters
  const double LAM_L = pars -> get_float("lambda_l");
  const double LAM_C = pars -> get_float("lambda_c"); 
  //hyp-smearing parameters
  const double ALPHA_1 = pars -> get_float("alpha_1");
  const double ALPHA_2 = pars -> get_float("alpha_2");
  const int ITER = pars -> get_int("iter");
  // set the chunks to calculate
  int size = MPI::COMM_WORLD.Get_size();
  int me = MPI::COMM_WORLD.Get_rank();
  int tstart = 0;
  int tend = 0;
  int todo = 0;
  int *tstarts = new int[size];
  int *tends = new int[size];
  int *todos = new int[size];
  if (me == 0) {
    // calculate chunks minimal chunk size
    int tmp = L0/size;
    // fill array
    for(int i = 0; i < size; ++i)
      todos[i] = tmp;
    // increase the first chunks until all
    // time slices are distributed
    for(int i = 0; i < L0%size; ++i)
      ++todos[i];
    // fill tstarts and tend;
    tstarts[0] = 0;
    tends[0] = todos[0] - 1;
    for(int i = 1; i < size; ++i) {
      tstarts[i] = tends[i-1] + 1;
      tends[i] = tstarts[i] + todos[i] - 1;
    }
    for(int i = 0; i < size; ++i) {
      std::cout << tstarts[i] << " " << tends[i] << std::endl;
    }
  }
  // scatter to all processes
  MPI::COMM_WORLD.Scatter(tstarts, 1, MPI_INT, &tstart, 1, MPI_INT, 0);
  MPI::COMM_WORLD.Scatter(tends, 1, MPI_INT, &tend, 1, MPI_INT, 0);
  MPI::COMM_WORLD.Scatter(todos, 1, MPI_INT, &todo, 1, MPI_INT, 0);

  printf("%d: size: %d, tstart %d, tend %d, todo %d\n", me, size, tstart, tend, todo);

  //lookup tables
  //N: # rows/columns, nev: desired # Eigenvalues, nconv: # converged EVs
  Tslice* slice = Tslice::getInstance();
  slice -> Tslice::init();
  PetscInt nev = NEV;
  PetscInt n;
  PetscInt nconv;
  PetscErrorCode ierr;
  //variables holding eigensystem
  PetscScalar eigr;
  PetscScalar eigi;
  std::vector<double> evals_accel;
  std::vector<double> phase;
  Vec xr;
  Vec xi;
  //Time Tracking
  PetscLogDouble v1;
  PetscLogDouble v2;
  PetscLogDouble elapsed_time;
  //Matrix to hold the entire eigensystem
  Eigen::MatrixXcd eigensystem(MAT_ENTRIES, nev);
  Eigen::MatrixXcd eigensystem_fix(MAT_ENTRIES, nev);
  std::complex<double> trc;

  //get number of configuration from last argument to main
  int config = atoi( argv[ (argc-1) ] );
  --argc;
  char conf_name [200];
  sprintf(conf_name, "%s/conf.%04d", GAUGE_FIELDS.c_str(), config );
  printf("%s\n", conf_name);
  printf("Using chebyshev parameters Lambda_c: %f and Lambda_l: %f\n", LAM_C, LAM_L);
  //Initialize memory for configuration
  double* configuration = new double[todo*V_TS];
  //double* configuration = new double[V_4_LIME];
  ierr = read_lime_gauge_field_doubleprec_timeslices(configuration, conf_name,
      L0, L1, L2, L3, tstart, tend);
  //__Initialize SLEPc__
  SlepcInitialize(&argc, &argv, (char*)0, NULL);
  std::cout << "initialized Slepc" << std::endl;
  //loop over timeslices of a configuration
  for (int ts = 0; ts < todo; ++ts) {
    std::cout << me << ": calculating time slice " << ts << std::endl;
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
    n = V3;//Tell Shell matrix about size of vectors
    std::cout << me << ": Try to create Shell Matrix..." << std::endl;
    ierr = MatCreateShell(PETSC_COMM_SELF,MAT_ENTRIES,MAT_ENTRIES,PETSC_DECIDE,
        PETSC_DECIDE,&n,&A);
      CHKERRQ(ierr);
    std::cout << me << ": done" << std::endl;
    std::cout << me << ": Try to set operations..." << std::endl;
    ierr = MatSetFromOptions(A);
      CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult_Laplacian2D);
      CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,
        (void(*)())MatMult_Laplacian2D);
      CHKERRQ(ierr);
    std::cout << me << ": accomplished" << std::endl;
    ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,
        (void(*)())MatGetDiagonal_Laplacian2D);
      CHKERRQ(ierr);
    std::cout << me << ": Matrix creation..." << std::endl;
      CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_HERMITIAN, PETSC_TRUE);
      CHKERRQ(ierr);
    ierr = MatSetUp(A);
      CHKERRQ(ierr);
    ierr = MatGetVecs(A,NULL,&xr);
      CHKERRQ(ierr);
    ierr = MatGetVecs(A,NULL,&xi);
      CHKERRQ(ierr);
      std::cout << me << ": successful" << std::endl;  

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

    //--------------------------------------------------------------------------//
    //                         Solve Ax = kx                                    //
    //--------------------------------------------------------------------------//

    ierr = PetscTime(&v1);
      CHKERRQ(ierr);
    std::cout << me << ": Start solving T_8(B) * x = y" << std::endl;
    ierr = EPSSolve(eps);
      CHKERRQ(ierr);
    ierr = PetscTime(&v2);
      CHKERRQ(ierr);
    elapsed_time = v2 - v1;

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
    ierr = PetscPrintf(PETSC_COMM_SELF,"Number of converged eigenvalues: %D\n",nconv);
    CHKERRQ(ierr);
    ierr=PetscPrintf(PETSC_COMM_SELF, "Convergence reason of solver: %D\n",reason);
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
    write_eig_sys_bin("eigenvectors", config, ts, nev, eigensystem_fix); 
    //recover spectrum of eigenvalues from acceleration
    std::vector<double> evals_save;
    recover_spectrum(nconv, evals_accel, evals_save);
    std::cout << evals_accel.at(0) << " " << evals_save.at(0) <<std::endl;
    write_eigenvalues_bin("eigenvalues", config, ts, nev, evals_save);
    write_eigenvalues_bin("phases", config, ts, nev, phase );
   
    //check trace of eigensystem
    trc = ( eigensystem.adjoint() * ( eigensystem ) ).trace();
    std::cout << "V.adj() * V = " << trc << std::endl;
    //Display user information on files
    printf("%d: %d phase fixed eigenvectors saved successfully \n", me, nconv);
    printf("%d: %d eigenvalues saved successfully \n", me, nconv);
    //__Clean up__
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
  } //end loop over timeslices

  ierr = SlepcFinalize();
  MPI::Finalize();
  return 0;
}


