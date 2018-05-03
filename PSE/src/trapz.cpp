#include "trapz.hpp"
#include "Init_Vec.hpp"
#include "set_Vec.hpp"
#include <iostream>

namespace PSE{
    int trapz(
            Vec &q,
            const PetscInt &n,
            PetscScalar &I,
            const PetscScalar &Deltax ///< Domain size (b-a) in \f[ I = \int_a^b f(x) dx\f]
            ){
        // set trapezoidal one row matrix
        Vec Trap;
        PetscErrorCode ierr;
        Init_Vec(Trap,n);
        PetscInt Istart, Iend;
        VecGetOwnershipRange(Trap,&Istart,&Iend);
        PetscScalar value= Deltax/((PetscScalar)n);
        PetscScalar valueBC = Deltax/(2.*(PetscScalar)n);
        
        for(int i=Istart; i<Iend; ++i){
            if(i==0 || i==n) {// set ends with half the value
                ierr = VecSetValue(Trap,i,valueBC,INSERT_VALUES);CHKERRQ(ierr);
            }
            else{ // set mid values
                ierr = VecSetValue(Trap,i,value,INSERT_VALUES);CHKERRQ(ierr);
            }
        }
        set_Vec(Trap);//assemble Trapezoidal matrix

        // I = Trap dot q
        VecDot(Trap,q,&I);

        VecDestroy(&Trap);

        return 0;
    }

}
