#include "base_flow.hpp"
int base_flow(PetscScalar *U,
        PetscScalar *Uy,
        PetscScalar Uyy,
        PetscScalar *y,
        int n,
        PetscBool output_full,
        PetscInt ){
    //PetscScalar U[n], Uy[n], Uyy;
    for(int i=0; i<n; i++){
        U[i] = 1. - y[i]*y[i];
        Uy[i]= -2.*y[i];
    }
    Uyy = -2.;

    return 0;
}
