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
            PSE &data
            ){
        // base flow
        if (data.flag_base_flow==PETSC_FALSE){
            base_flow(data.U,data.Uy,data.y,data.ny);
            data.flag_base_flow=PETSC_TRUE;
        }
        // set D operators
        if (data.flag_set_D==PETSC_FALSE){
            set_D(data.y,data.ny,data.Dy,data.order,1);
            //printMatASCII(Dy,"Dy.txt");
            set_D(data.z,data.nz,data.Dz,data.order,1,PETSC_TRUE);
            set_D(data.z,data.nz,data.Dzz,data.order,2,PETSC_TRUE);
            set_D(data.y,data.ny,data.Dyy,data.order,2);
            data.flag_set_D=PETSC_TRUE;
        }
        // zero A and B matrix
        
        for (PetscInt i=0; i<data.dim; ++i){
            MatZeroRows(data.A,
        }

        //Vec b;
        PetscErrorCode ierr;

        for (PetscInt i=0; i<data.nz; i++){
            PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_before_assemble set zi=%d",i);
            set_A_and_B_zi(data.y,data.ny,data.z,data.nz,data.A,data.B,data.Re,data.rho,data.alpha,data.m,data.omega,data.Dy,data.Dyy,data.Dz,data.Dzz,data.I,data.U,data.Uy,i,data.order,data.reduce_wall_order);
        }
        // assemble A,B
        PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_before_assemble");
        set_Mat(data.A);
        PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_in_assemble");
        set_Mat(data.B);
        PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_after_assemble");


        //ierr = MatDestroy(&B);CHKERRQ(ierr);

        return 0;
    }
}
