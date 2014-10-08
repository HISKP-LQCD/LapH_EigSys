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
#include <string>
#include <Eigen/Core>
#include "slepceps.h"
#include "config_utils.h"
#include "eigensystem.h"
#include "io.h"
#include "par_io.h"
#include "shell_matop.h"
#include "variables.h"
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

  Mat A;              
  EPS eps;             
  EPSType type;
  EPSConvergedReason reason;
  //Handling infile
  IO* pars = IO::getInstance();
  pars -> set_values("parameters.txt");
  pars -> print_summary();
  //in and outpaths
  std::string GAUGE_FIELDS = pars -> get_path("config_path");
  //lattice layout from infile 
  int L0 = pars -> get_int("LT");
  int L1 = pars -> get_int("LX");
  int L2 = pars -> get_int("LY");
  int L3 = pars -> get_int("LZ");
  const int V3 = pars -> get_int("V3");
  int V_TS = pars -> get_int("V_TS");
  //calculation parameters from infile
  int NEV = pars -> get_int("NEV");
  int V_4_LIME = pars -> get_int("V4_LIME");
  const int MAT_ENTRIES = pars -> get_int("MAT_ENTRIES");
  //chebyshev parameters
  int LAM_L = pars -> get_int("LAM_L");
  int LAM_C = pars -> get_int("LAM_C"); 
  //hyp-smearing parameters
  double ALPHA_1 = pars -> get_float("alpha_1");
  double ALPHA_2 = pars -> get_float("alpha_2");
  int ITER = pars -> get_int("iter");

  //lookup tables
  int up_3d[V3][3],down_3d[V3][3];
  //Eigen Array
  Eigen::Matrix3cd **eigen_timeslice = new Eigen::Matrix3cd *[V3];
  //N: # rows/columns, nev: desired # Eigenvalues, nconv: # converged EVs
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

  //Allocate Eigen Array to hold timeslice
  for ( auto i = 0; i < V3; ++i ) {
    eigen_timeslice[i] = new Eigen::Matrix3cd[3];
    for (auto dir = 0; dir < 3; ++dir) {
      eigen_timeslice[i][dir] = Eigen::Matrix3cd::Identity();
    }
  }
  //Initialize lookup-tables
  hopping3d(up_3d, down_3d);
  //set up output
  //get number of configuration from last argument to main
  int config;
  config = atoi( argv[ (argc-1) ] );
  --argc;
  char conf_name [200];
  sprintf(conf_name, "%s/conf.%04d", GAUGE_FIELDS.c_str(), config );
  printf("%s\n", conf_name);
  printf("Using chebyshev parameters Lambda_c: %f and Lambda_l: %f\n", LAM_C, LAM_L);
  //Initialize memory for configuration
  double* configuration = new double[V_4_LIME];
  ierr = read_lime_gauge_field_doubleprec_timeslices(configuration, conf_name,
      L0, L1, L2, L3, 0, L0);
  //__Initalize SLEPc__
  SlepcInitialize(&argc, &argv, (char*)0, NULL);
  //loop over timeslices of a configuration
  for (int ts = 0; ts < L0; ++ts) {

    //--------------------------------------------------------------------------//
    //                              Data input                                  //
    //--------------------------------------------------------------------------//

    //for every timeslice open new file binary write
    //evectors = fopen(ts_name,"wb");
    //Time Slice of Configuration
    double* timeslice = configuration + (ts*V_TS);

    //Write Timeslice in Eigen Array                                                  
    map_timeslice_to_eigen(eigen_timeslice, timeslice);

    //Apply Smearing algorithm to timeslice ts
    smearing_hyp(eigen_timeslice,ALPHA_1, ALPHA_2, ITER);
    //__Define Action of Laplacian in Color and spatial space on vector
    n = V3;//Tell Shell matrix about size of vectors
    ierr = MatCreateShell(PETSC_COMM_WORLD,MAT_ENTRIES,MAT_ENTRIES,PETSC_DECIDE,
        PETSC_DECIDE,&n,&A);
      CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
      CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult_Laplacian2D);
      CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,
        (void(*)())MatMult_Laplacian2D);
      CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,
        (void(*)())MatGetDiagonal_Laplacian2D);
      CHKERRQ(ierr);
    std::cout << "Matrix creation: ";
      CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_HERMITIAN, PETSC_TRUE);
      CHKERRQ(ierr);
    ierr = MatSetUp(A);
      CHKERRQ(ierr);
    ierr = MatGetVecs(A,NULL,&xr);
      CHKERRQ(ierr);
    ierr = MatGetVecs(A,NULL,&xi);
      CHKERRQ(ierr);
      std::cout << "successful" << std::endl;  
    std::cout << "Start solving T_8(B) * x = y" << std::endl;

    //--------------------------------------------------------------------------//
    //                 Context creation & Options setting                       //
    //--------------------------------------------------------------------------//

    //Create eigensolver context
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
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
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of converged eigenvalues: %D\n",nconv);
    CHKERRQ(ierr);
    ierr=PetscPrintf(PETSC_COMM_WORLD, "Convergence reason of solver: %D\n",reason);
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
    
      //save nconv eigenvectors and eigenvalues to one file each
      //fwrite(ptr_evec, sizeof(ptr_evec[0]), MAT_ENTRIES, evectors);
      //write each eigenvector to eigensystem for check
      eigensystem.col(i) = Eigen::Map<Eigen::VectorXcd, 0 >(ptr_evec, MAT_ENTRIES);
      //fwrite(ptr_eval ,  sizeof(double) , 1, evalues );
      evals_accel.at(i) = eval;
      ierr = VecRestoreArray(xr, &ptr_evec);
      CHKERRQ(ierr);
    }
    //fix phase to 0 in first entry of each eigenvector
    fix_phase(eigensystem, eigensystem_fix, phase);
    write_eig_sys_bin("eigenvectors", config, ts, nev, eigensystem_fix); 
    //recover spectrum of eigenvalues from acceleration
    std::vector<double> evals_save;
    recover_spectrum(nconv, evals_accel, evals_save);
    write_eigenvalues_bin("eigenvalues", config, ts, nev, evals_save);
    write_eigenvalues_bin("phases", config, ts, nev, phase );
    /*
    eigenvalues.write(reinterpret_cast<char*>(&evals_accel[0]), evals_accel.size()*sizeof(double));
    */
    //check trace of eigensystem
    trc = ( eigensystem.adjoint() * ( eigensystem ) ).trace();
    std::cout << "V.adj() * V = " << trc << std::endl;
    //Display user information on files
    printf("%d eigenvectors saved successfully \n", nconv);
    printf("%d eigenvalues saved successfully \n", nconv);
    //__Clean up__
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
  }//end loop over timeslices

  ierr = SlepcFinalize();
  //delete configuration;
  for (auto j = 0; j < V3; ++j) delete[] eigen_timeslice[j];
  delete eigen_timeslice;
  return 0;
}


