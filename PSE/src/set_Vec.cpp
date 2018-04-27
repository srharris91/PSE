#include "set_Vec.hpp"
#include "Init_Vec.hpp"
#include <iostream>

namespace PSE{
    PetscInt set_Vec(
            const PetscScalar *bin,
            const PetscInt &n,
            Vec &b,
            const PetscBool &parallel
            ){
        PetscInt Istart, Iend, Ii;
        PetscErrorCode ierr;
        // set Vec in parallel
        if (parallel){
            ierr = VecGetOwnershipRange(b,&Istart,&Iend);CHKERRQ(ierr);
            for (Ii=Istart; Ii<Iend; Ii++){ // in parallel
                ierr = VecSetValues(b,1,&Ii,&bin[Ii],INSERT_VALUES);CHKERRQ(ierr);
            }

            // broadcast
            ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
        }
        else{
            for (Ii=0; Ii<n; Ii++){
                ierr = VecSetValues(b,1,&Ii,&bin[Ii],INSERT_VALUES);CHKERRQ(ierr);
            }
        }

        return 0;
    }
    PetscInt set_Vec(
            Vec &b
            ){
        PetscErrorCode ierr;
        ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
        return 0;
    }
}

