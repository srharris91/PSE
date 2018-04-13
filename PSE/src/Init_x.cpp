#include "Init_x.hpp"
#include "print.hpp"
#include <iostream>
namespace PSE{
    int Init_x(
            Vec &x, 
            PetscInt &n
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
