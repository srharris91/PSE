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
#include "class.hpp"

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
        MatDestroy(&data.A);
        MatDestroy(&data.B);
        Init_Mat(data.A,data.dim);
        Init_Mat(data.B,data.dim);
        /*
        PetscInt *temp_row = new PetscInt[data.dim];
        for (PetscInt i=0; i<data.dim; ++i) temp_row[i] = i; // set row numbers
        MatZeroRows(data.A,data.dim,temp_row,0,0,0);
        delete[] temp_row;
        */

        //Vec b;
        PetscErrorCode ierr;

        for (PetscInt i=0; i<data.nz; i++){
            set_A_and_B_zi(data.y,data.ny,data.z,data.nz,data.A,data.B,data.Re,data.rho,data.alpha,data.m,data.omega,data.Dy,data.Dyy,data.Dz,data.Dzz,data.I,data.U,data.Uy,i,data.order,data.reduce_wall_order);
        }
        // assemble A,B
        set_Mat(data.A);
        set_Mat(data.B);


        //ierr = MatDestroy(&B);CHKERRQ(ierr);

        return 0;
    }
}
