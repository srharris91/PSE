#include "base_flow.hpp"
#include "Init_Mat.hpp"
#include "set_Mat.hpp"

namespace PSE
{
    int base_flow(
            PetscScalar U[],
            PetscScalar Uy[],
            const PetscScalar y[],
            const PetscInt &n
            ){
        // set U,Uy
        for(int i=0; i<n; i++){
            U[i] = 1. - y[i]*y[i];
            Uy[i]= -2.*y[i];
        }
        return 0;
    }

    int base_flow(
            Mat &U,
            Mat &Uy,
            const PetscScalar y[],
            const PetscInt &n
            ){
        // init U,Uy
        Init_Mat(U,n);
        Init_Mat(Uy,n);
        // set U,Uy
        for(int i=0; i<n; i++){
            set_Mat(1. - y[i]*y[i],i,i,U);
            set_Mat(-2.*y[i],i,i,Uy);
        }
        // assemble U,Uy
        set_Mat(U);
        set_Mat(Uy);

        return 0;
    }
}
