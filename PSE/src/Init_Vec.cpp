#include "Init_Vec.hpp"
#include "print.hpp"
namespace PSE{
    int Init_Vec(
            Vec &x, 
            const PetscInt &n
            ){
        PetscErrorCode ierr;
        ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
        ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
        ierr = VecSetFromOptions(x);CHKERRQ(ierr);
        //ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        //ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
        return 0;
    }
}
