#include <stdexcept>
#include <string>
#include <cmath>
#include "set_A_and_B.hpp"
#include "set_D.hpp"
#include "get_D_Coeffs.hpp"
#include "Init_Vec.hpp"
#include "Init_Mat.hpp"
#include "print.hpp"
#include "set_Mat.hpp"
#include "set_Vec.hpp"
#include <iostream>

namespace PSE
{
    PetscInt set_A_and_B(
            const PetscScalar y[],
            const PetscInt &ny,
            const PetscScalar z[],
            const PetscInt &nz,
            Mat &A,         
            Mat &B,         
            const PetscInt &order,
            const PetscBool &reduce_wall_order
            ){
        PetscErrorCode ierr;
        PetscInt dim=ny*nz*4; // dimension of A and B square matrices
        // initialize D operators
        Mat Dy,Dyy,Dz,Dzz;
        set_D(y,ny,Dy,order,1);
        set_D(z,nz,Dz,order,1);
        set_D(z,nz,Dzz,order,2);
        set_D(y,ny,Dyy,order,2);
        
        Init_Mat(A,dim);
        Init_Mat(B,dim);

        // set u-mom equation
        set_Mat(Dy,ny,B,dim,ny,ny);
        set_Mat(Dy,ny,B,dim);
        set_Mat(Dz,nz,ny,0,A,dim);

        // assemble
        set_Mat(A);
        set_Mat(B);
        // destroy
        ierr = MatDestroy(&Dy);CHKERRQ(ierr);
        ierr = MatDestroy(&Dyy);CHKERRQ(ierr);
        ierr = MatDestroy(&Dz);CHKERRQ(ierr);
        ierr = MatDestroy(&Dzz);CHKERRQ(ierr);

        return 0;
    }
}
