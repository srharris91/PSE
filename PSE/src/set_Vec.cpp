#include "set_Vec.hpp"
#include "Init_Vec.hpp"

namespace PSE{
    PetscInt set_Vec(
            const PetscScalar *bin,
            const PetscInt &n,
            Vec &b
            ){
        PetscInt Istart, Iend, Ii;
        PetscErrorCode ierr;
        // set Vec
        ierr = VecGetOwnershipRange(b,&Istart,&Iend);CHKERRQ(ierr);
        for (Ii=Istart; Ii<Iend; Ii++){ // in parallel
            ierr = VecSetValues(b,1,&Ii,&bin[Ii],INSERT_VALUES);CHKERRQ(ierr);
        }

        // broadcast
        ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

        return 0;
        }
    }

