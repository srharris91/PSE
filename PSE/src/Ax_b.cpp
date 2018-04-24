#include "Init_Mat.hpp"
#include "Init_Vec.hpp"
#include "Ax_b.hpp"
#include "print.hpp"
#include "set_Mat.hpp"
#include "set_Vec.hpp"

namespace PSE
{
    PetscInt Ax_b(
            const PetscScalar* const* Ain,
            Vec &x,
            const PetscScalar bin[],
            const PetscInt &n 
            ){
        Vec            b;      /* approx solution, RHS, exact solution */
        Mat            A;            /* linear system matrix */
        KSP            ksp;         /* linear solver context */
        PetscErrorCode ierr;
        
        // set A matrix
        set_Mat(Ain,n,A);
        set_Mat(A);

        // Create parallel vectors.
        ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

        // set right-hand-side vector.
        set_Vec(bin,n,b);

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

        // Free work space.  All PETSc objects should be destroyed when they
        ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
        // are no longer needed.
        ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
        //ierr = VecDestroy(&x);CHKERRQ(ierr);
        ierr = VecDestroy(&b);CHKERRQ(ierr); 
        ierr = MatDestroy(&A);CHKERRQ(ierr);
        //ierr = PetscFinalize();
        return 0;
    }
    PetscInt Ax_b(
            const Mat &A,
            Vec &x,
            const Vec &b,
            const PetscInt &n 
            ){
        KSP            ksp;         /* linear solver context */
        PetscErrorCode ierr;
        
        // Create parallel vectors x from b
        Init_Vec(x,n);
        //ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
        //ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

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

        // output iterations
        PetscInt its;
        ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"ksp iterations %D\n",its);CHKERRQ(ierr);
        // Free work space.  All PETSc objects should be destroyed when they
        ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
        // are no longer needed.
        ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,"KSP Solved\n");
        //ierr = VecDestroy(&x);CHKERRQ(ierr);
        //ierr = VecDestroy(&b);CHKERRQ(ierr); 
        //ierr = MatDestroy(&A);CHKERRQ(ierr);
        //ierr = PetscFinalize();
        return 0;
    }
}
