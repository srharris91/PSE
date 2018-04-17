#include "set_Mat.hpp"
#include "Init_Mat.hpp"

namespace PSE{
    PetscInt set_Mat(
            const PetscScalar* const* Ain,
            const PetscInt &n,
            Mat &A
            ){
        PetscInt Istart, Iend, Ii, col[n],i;
        for(i=0;i<n;i++) col[i]=i;
        PetscErrorCode ierr;

        // initialize A matrix of right size
        Init_Mat(A,n);

        // Currently, all PETSc parallel matrix formats are partitioned by
        // contiguous chunks of rows across the processors.  Determine which
        // rows of the matrix are locally owned.
        ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

        // Set matrix elements (in parallel)

        for (Ii=Istart; Ii<Iend; Ii++){ // in parallel all of the values
            ierr = MatSetValues(A,1,&Ii,n,col,Ain[Ii],INSERT_VALUES);CHKERRQ(ierr);
        }


        // Assemble matrix, using the 2-step process:
        // MatAssemblyBegin(), MatAssemblyEnd()
        // Computations can be done while messages are in transition
        // by placing code between these two statements.
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

        return 0;
    }
    PetscInt set_Mat(
            const PetscScalar &diag,
            const PetscInt &n,
            Mat &A,
            const PetscInt &k,
            const PetscBool &parallel
            ){
        PetscErrorCode ierr;
        PetscInt Istart, Iend, Ii;

        // initialize A matrix of right size
        //Init_Mat(A,n);

        // Currently, all PETSc parallel matrix formats are partitioned by
        // contiguous chunks of rows across the processors.  Determine which
        // rows of the matrix are locally owned.
        ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

        // Set matrix elements (in parallel)

        if (parallel){
            for (Ii=Istart; Ii<Iend; Ii++){ // in parallel all of the values
                if (Ii+k >= 0 && Ii+k <= n-1){
                    ierr = MatSetValue(A,Ii,(Ii + k),diag,INSERT_VALUES);CHKERRQ(ierr);
                }
            }
        }
        else{
            for (Ii=0; Ii<n; Ii++){ // in parallel all of the values
                if (Ii+k >= 0 && Ii+k <= n-1){
                    ierr = MatSetValue(A,Ii,(Ii + k),diag,INSERT_VALUES);CHKERRQ(ierr);
                }
            }
        }


        return 0;
    }
    PetscInt set_Mat(
            Mat &A  //< Mat to assemble
            ){
        // Assemble matrix, using the 2-step process:
        // MatAssemblyBegin(), MatAssemblyEnd()
        // Computations can be done while messages are in transition
        // by placing code between these two statements.

        PetscErrorCode ierr;
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        return 0;
    }
}

