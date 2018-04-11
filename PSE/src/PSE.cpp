
static char help[] = "Solves a linear system in parallel with KSP.\n\n";

/*T
   Concepts: KSP^solving the PSE equations
   Concepts: complex numbers;
   Concepts: Parabolic Stability Equations for Transition modeling
   Processors: n
T*/

/*
   Description: Solves a complex linear system in parallel with KSP.

   The model problem:
      Solve PSE equation in the channel flow
      Dirichlet b.c.'s on all sides
      Use the 3-D, finite difference stencil.

   Compiling the code:
      This code uses the complex numbers version of PETSc, so configure
      must be run to enable this
*/

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include <petscksp.h>
#include <stdio.h>

PetscScalar* base_flow(PetscScalar *y, PetscScalar n, PetscBool output_full=PETSC_FALSE, PetscInt Dim=2){
    PetscScalar U[int(n)], Uy[int(n)], Uyy;
    for(int i=0; i<n; i++){
        U[i] = 1. - y[i]*y[i];
        Uy[i]= -2.*y[i];
    }
    Uyy = -2.;

    return 0;
}
int Read_q(PetscScalar *RHS_True,int n){
    PetscErrorCode ierr;
    char buff[]="tofile.dat";
    FILE *latfile;
    latfile=fopen(buff,"r");
    fread(RHS_True,sizeof(double),n,latfile);
    fclose(latfile);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nReading in RHS_True\n");CHKERRQ(ierr);
    for(int i=0; i<n; i++){
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  RHS_True[%D] = %g \n",i,(double)RHS_True[i]);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"done Reading in RHS_True\n\n\n");CHKERRQ(ierr);
    return 0;

}

int main(int argc,char **args)
{
  Vec            x,b,u;      /* approx solution, RHS, exact solution */
  Mat            A;            /* linear system matrix */
  KSP            ksp;         /* linear solver context */
  PetscReal      norm;         /* norm of solution error */
  PetscInt       dim,i,Istart,Iend,n = 6,its,col[]={0,1,2,3,4,5};//,j
  PetscErrorCode ierr;
  PetscScalar    none = -1.0,*xa,v[n],RHS_True[n];
  PetscBool      flg = PETSC_FALSE;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  dim  = n*n;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,dim,dim);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  /*
     Set matrix elements 
  */
  
  v[0]=1;   v[1]=1; v[2]=-2;v[3]=1; v[4]=3; v[5]=-4.;  i=0; ierr = MatSetValues(A,1,&i,n,col,v,INSERT_VALUES);CHKERRQ(ierr);
  v[0]=2;   v[1]=-1;v[2]=1; v[3]=2; v[4]=1; v[5]=-3.;  i=1; ierr = MatSetValues(A,1,&i,n,col,v,INSERT_VALUES);CHKERRQ(ierr);
  v[0]=1;   v[1]=3; v[2]=-3;v[3]=-1;v[4]=2; v[5]=1. ;  i=2; ierr = MatSetValues(A,1,&i,n,col,v,INSERT_VALUES);CHKERRQ(ierr);
  v[0]=5;   v[1]=2; v[2]=-1;v[3]=-1;v[4]=2; v[5]=1. ;  i=3; ierr = MatSetValues(A,1,&i,n,col,v,INSERT_VALUES);CHKERRQ(ierr);
  v[0]=-3;  v[1]=-1;v[2]=2; v[3]=3; v[4]=1; v[5]=3. ;  i=4; ierr = MatSetValues(A,1,&i,n,col,v,INSERT_VALUES);CHKERRQ(ierr);
  v[0]=4;   v[1]=3; v[2]=1; v[3]=-6;v[4]=-3;v[5]=-2.;  i=5; ierr = MatSetValues(A,1,&i,n,col,v,INSERT_VALUES);CHKERRQ(ierr);


  /*
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd()
     Computations can be done while messages are in transition
     by placing code between these two statements.
  */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /*
     Create parallel vectors.
      - When using VecCreate(), VecSetSizes() and VecSetFromOptions(),
      we specify only the vector's global
        dimension; the parallel partitioning is determined at runtime.
      - Note: We form 1 vector from scratch and then duplicate as needed.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,dim);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

  /*
     Set exact solution; then compute right-hand-side vector.
  */

  //ierr = VecSet(u,pfive);CHKERRQ(ierr);
  Read_q(RHS_True,n);
  ierr = VecSetValues(u,n,col,RHS_True,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatMult(A,u,b);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
      Print the first 3 entries of x; this demonstrates extraction of the
      real and imaginary components of the complex vector, x.
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-print_x6",&flg,NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nThe first six entries of x are:\n");CHKERRQ(ierr);
    for (i=0; i<6; i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  x[%D] = %g + %g i\n",i,(double)PetscRealPart(xa[i]),(double)PetscImaginaryPart(xa[i]));CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n");CHKERRQ(ierr);
    ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);
  }

  /*
     Check the error
  */
  ierr = VecAXPY(x,none,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  if (norm < 1.e-12) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error < 1.e-12 iterations %D\n",its);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %D\n",(double)norm,its);CHKERRQ(ierr);
  }

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr); ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
