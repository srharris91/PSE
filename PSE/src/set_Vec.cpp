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
    PetscInt set_Vec(
            const PetscScalar &bin,
            const PetscInt &n,
            Vec &b
            ){
        PetscErrorCode ierr;
        //set vec
        ierr = VecSetValues(b,1,&n,&bin,INSERT_VALUES);CHKERRQ(ierr);

        return 0;
    }
    PetscInt set_Vec(
            const Vec &inVec,
            const PetscInt &low,
            const PetscInt &hi,
            Vec &subVec
            ){
        // init subvec
        //Init_Vec(subVec,hi-low);
        // range for inVec
        PetscInt low_local,high_local;
        VecGetOwnershipRange(inVec,&low_local,&high_local);
        PetscScalar to_insert;
        for(int i=low_local; i<high_local; ++i) {
            if (i>=low && i<hi) {
                //std::cout<<"i = "<<i<<std::endl;
                VecGetValues(inVec,1,&i,&to_insert);
                //std::cout<<"to_insert = "<<to_insert<<std::endl;
                PetscInt loc=i-low;
                set_Vec(to_insert,loc,subVec);
            }

        }
        set_Vec(subVec);// assemble

        return 0;
    }
    
}

