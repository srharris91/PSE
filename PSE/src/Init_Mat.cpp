#include "Init_Mat.hpp"
namespace PSE{
    int Init_Mat(
            Mat &A, 
            const PetscInt &n
            ){
        PetscErrorCode ierr;
        ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
        ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
        ierr = MatSetFromOptions(A);CHKERRQ(ierr);
        ierr = MatSetUp(A);CHKERRQ(ierr);
        return 0;
    }
    int Init_Mat(
            Mat &A, 
            const PetscInt &m,
            const PetscInt &n
            ){
        PetscErrorCode ierr;
        ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
        ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);CHKERRQ(ierr);
        ierr = MatSetFromOptions(A);CHKERRQ(ierr);
        ierr = MatSetUp(A);CHKERRQ(ierr);
        return 0;
    }
}
