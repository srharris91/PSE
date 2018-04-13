#include "Ax_b.hpp"
#include <iostream>

namespace PSE
{
    PetscErrorCode Ax_b(
            PetscScalar **Ain,
            Vec &x,
            PetscScalar bin[],
            PetscInt n 
            ){
        Vec            b;      /* approx solution, RHS, exact solution */
        Mat            A;            /* linear system matrix */
        KSP            ksp;         /* linear solver context */
        PetscInt       dim,i,Ii,Istart,Iend,col[n];//,j
        for(i=0;i<n;i++) col[i]=i;
        PetscErrorCode ierr;
        //PetscScalar    *xa;
        //PetscScalar (*Ain)[n] = (PetscScalar (*)[n]) Ain1;

        //ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
        //ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
        dim  = n;

        // Compute the matrix and right-hand-side vector that define
        // the linear system, Ax = b.
        // Create parallel matrix, specifying only its global dimensions.
        // When using MatCreate(), the matrix format can be specified at
        // runtime. Also, the parallel partitioning of the matrix is
        // determined by PETSc at runtime.
        ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
        ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,dim,dim);CHKERRQ(ierr);
        ierr = MatSetFromOptions(A);CHKERRQ(ierr);
        ierr = MatSetUp(A);CHKERRQ(ierr);

        // Currently, all PETSc parallel matrix formats are partitioned by
        // contiguous chunks of rows across the processors.  Determine which
        // rows of the matrix are locally owned.
        ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

        // Set matrix elements (in parallel)

        for (Ii=Istart; Ii<Iend; Ii++){ // in parallel
            //std::cout<<"Ii = "<<Ii<<" of start "<<Istart<<" to end "<<Iend<<std::endl;
        //for(i=0; i<n; i++){  // in serial

            ierr = MatSetValues(A,1,&Ii,n,col,Ain[Ii],INSERT_VALUES);CHKERRQ(ierr);
        }


        // Assemble matrix, using the 2-step process:
        // MatAssemblyBegin(), MatAssemblyEnd()
        // Computations can be done while messages are in transition
        // by placing code between these two statements.
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

        // Create parallel vectors.
        // - When using VecCreate(), VecSetSizes() and VecSetFromOptions(),
        // we specify only the vector's global
        // dimension; the parallel partitioning is determined at runtime.
        // - Note: We form 1 vector from scratch and then duplicate as needed.
        //ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        //ierr = VecSetSizes(b,PETSC_DECIDE,dim);CHKERRQ(ierr);
        //ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

        // set right-hand-side vector.
        ierr = VecSetValues(b,n,col,bin,INSERT_VALUES);CHKERRQ(ierr);

        ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

        //ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        //ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

        // Create the linear solver and set various options
        // Create linear solver context
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        // Set operators. Here the matrix that defines the linear system
        // also serves as the preconditioning matrix.
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        // Set runtime options, e.g.,
        // -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

        // Solve the linear system
        ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

        // Check solution and clean up
        // Print the first 3 entries of x; this demonstrates extraction of the
        // real and imaginary components of the complex vector, x.
        //ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nThe first six entries of x are:\n");CHKERRQ(ierr);
        //for (i=0; i<6; i++) {
            //ierr = PetscPrintf(PETSC_COMM_WORLD,"  x[%D] = %g + %g i\n",i,(double)PetscRealPart(xa[i]),(double)PetscImaginaryPart(xa[i]));CHKERRQ(ierr);
            //xout[i] = xa[i]; // copy to xout
        //}
        //ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n");CHKERRQ(ierr);
        //ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);

        // Free work space.  All PETSc objects should be destroyed when they
        // are no longer needed.
        ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
        //ierr = VecDestroy(&x);CHKERRQ(ierr);
        ierr = VecDestroy(&b);CHKERRQ(ierr); 
        ierr = MatDestroy(&A);CHKERRQ(ierr);
        //ierr = PetscFinalize();
        return ierr;
    }
}
