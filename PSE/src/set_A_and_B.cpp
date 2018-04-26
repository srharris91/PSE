#include <stdexcept>
#include <string>
#include <cmath>
#include "set_A_and_B_zi.hpp"
#include "set_A_and_B.hpp"
#include "set_D.hpp"
#include "get_D_Coeffs.hpp"
#include "Init_Vec.hpp"
#include "Init_Mat.hpp"
#include "print.hpp"
#include "set_Mat.hpp"
#include "set_Vec.hpp"
#include "base_flow.hpp"

namespace PSE
{
    PetscInt set_A_and_B(
            const PetscScalar y[],
            const PetscInt &ny,
            const PetscScalar z[],
            const PetscInt &nz,
            Mat &A,         
            Mat &B,         
            const PetscScalar &Re,
            const PetscScalar &rho,
            const PetscScalar &alpha,
            const PetscScalar &m,
            const PetscScalar &omega,
            const PetscInt &order,
            const PetscBool &reduce_wall_order
            ){
        // calc base flow and derivatives
        // initialize D operators
        Mat Dy,Dyy,Dz,Dzz,I,U,Uy;
        // set I
        Init_Mat(I,ny);
        set_Mat(1.,ny,I); // set diag of 1
        set_Mat(I); // assemble
        // set U,Uy
        base_flow(U,Uy,y,ny);
        // set D operators
        set_D(y,ny,Dy,order,1);
        set_D(z,nz,Dz,order,1,PETSC_TRUE);
        set_D(z,nz,Dzz,order,2,PETSC_TRUE);
        set_D(y,ny,Dyy,order,2);

        //Vec b;
        PetscErrorCode ierr;
        PetscInt dim=ny*nz*4;

        Init_Mat(A,dim);
        Init_Mat(B,dim);
        //Init_Vec(b,dim);
        for (PetscInt i=0; i<nz; i++){
            set_A_and_B_zi(y,ny,z,nz,A,B,Re,rho,alpha,m,omega,Dy,Dyy,Dz,Dzz,I,U,Uy,i,order,reduce_wall_order);
        }
        // assemble A,B
        set_Mat(A);
        set_Mat(B);


        //ierr = MatDestroy(&B);CHKERRQ(ierr);
        ierr = MatDestroy(&Dy);CHKERRQ(ierr);
        ierr = MatDestroy(&Dyy);CHKERRQ(ierr);
        ierr = MatDestroy(&Dz);CHKERRQ(ierr);
        ierr = MatDestroy(&Dzz);CHKERRQ(ierr);
        ierr = MatDestroy(&I);CHKERRQ(ierr);
        ierr = MatDestroy(&U);CHKERRQ(ierr);
        ierr = MatDestroy(&Uy);CHKERRQ(ierr);

        return 0;
    }
}
