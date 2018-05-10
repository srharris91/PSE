#include "trapz.hpp"
#include "Init_Vec.hpp"
#include "set_Vec.hpp"
#include <iostream>

namespace PSE{
    int trapz(
            const Vec &q,
            const PetscInt &n,
            PetscScalar &I,
            const PetscScalar &Deltay 
            ){
        // set trapezoidal column (one row matrix)
        Vec Trap;
        PetscErrorCode ierr;
        Init_Vec(Trap,n);
        PetscInt Istart, Iend;
        VecGetOwnershipRange(Trap,&Istart,&Iend);
        PetscScalar value= Deltay/((PetscScalar)n);
        PetscScalar valueBC = Deltay/(2.*(PetscScalar)n);
        
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
    int trapz(
            const Vec &q,
            const PetscInt &ny,
            const PetscInt &nz,
            PetscScalar &I,
            const PetscScalar &Deltay,
            const PetscScalar &Deltaz
            ){
        // set trapezoidal column (one row matrix)
        Vec Trap;
        PetscErrorCode ierr;
        Init_Vec(Trap,ny);
        PetscInt Istart, Iend;
        VecGetOwnershipRange(Trap,&Istart,&Iend);
        PetscScalar value= Deltay/((PetscScalar)ny);
        PetscScalar valueBC = Deltay/(2.*(PetscScalar)ny);
        
        for(int i=Istart; i<Iend; ++i){
            if(i==0 || i==ny) {// set ends with half the value
                ierr = VecSetValue(Trap,i,valueBC,INSERT_VALUES);CHKERRQ(ierr);
            }
            else{ // set mid values
                ierr = VecSetValue(Trap,i,value,INSERT_VALUES);CHKERRQ(ierr);
            }
        }
        set_Vec(Trap);//assemble Trapezoidal matrix for y integration

        // I = Trap dot q
        Vec qsub,qsub_z;
        Init_Vec(qsub,ny);
        Init_Vec(qsub_z,nz);
        for (int zi=0; zi<nz; ++zi){ // for each z plane
            PetscScalar yi_uvwp = 0;
            VecZeroEntries(qsub);
            PetscScalar trap_value;
            for (int vari=0; vari<4; ++vari){       // for each variable u,v,w,P
                PetscInt row_eq = (4*zi + vari)*ny;
                set_Vec(q,row_eq,row_eq+ny,qsub);   // extract a single variable, all y values on z-plane
                VecDot(Trap,qsub,&trap_value);      // integrate in y
                yi_uvwp += trap_value;
            }
            set_Vec(yi_uvwp,zi,qsub_z); // set integral value into qsub_z vector
        }
        set_Vec(qsub_z);// assemble qsub_z
        trapz(qsub_z,nz,I,Deltaz); // integrate all values in z direction


        VecDestroy(&qsub);
        VecDestroy(&qsub_z);
        VecDestroy(&Trap);

        return 0;
    }

}
