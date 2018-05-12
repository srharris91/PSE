#include <stdexcept>
#include <string>
#include <cmath>
#include "set_A_and_B_zi.hpp"
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
    PetscInt set_A_and_B_zi(
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
            const Mat &Dy,
            const Mat &Dyy,
            const Mat &Dz,
            const Mat &Dzz,
            const Mat &I,
            const Mat &U,
            const Mat &Uy,
            const PetscInt &zi,
            const PetscInt &order,
            const PetscBool &reduce_wall_order
            ){
        PetscInt dim=ny*nz*4; // dimension of A and B square matrices
        PetscScalar i = 1.*PETSC_i; // imaginary number

        PetscInt row_eq;// row for the equation to be added (momentum u,v,w and continuity)
        PetscInt col_u=(4*zi+0)*ny, // set columns for variables in matrix
                 col_v=(4*zi+1)*ny,
                 col_w=(4*zi+2)*ny,
                 col_P=(4*zi+3)*ny;
        // set u-mom equation for A
        row_eq=(4*zi+0)*ny;
        PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_zi before terms");
        // u-terms
        set_Mat(-i*m*omega,I,ny,A,dim,row_eq,col_u);//-imw
        set_Mat(i*alpha,U,ny,A,dim,row_eq,col_u);//ialpha*U
        set_Mat((1./Re) * (pow(alpha,2)),I,ny,A,dim,row_eq,col_u);//Re^-1 * alpha^2
        PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_zi before most terms before Dyy");
        set_Mat(-1./Re,Dyy,ny,A,dim,row_eq,col_u);// -Re^-1 Dyy
        PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_zi before most terms after Dyy");
        set_Mat(-1./Re,Dzz,nz,ny,zi,A,dim,row_eq,col_u-4*zi*ny); // -Re^-1 Dzz
        PetscPrintf(PETSC_COMM_WORLD,"made it to here A_and_B_zi before most terms after Dzz");
        // v-terms
        set_Mat(1.,Uy,ny,A,dim,row_eq,col_v);//Uy
        //w-terms
        //P-terms
        set_Mat((1./rho) * i*alpha,I,ny,A,dim,row_eq,col_P);//imw
        // set u-mom equation for B
        // u-terms
        set_Mat(1.,U,ny,B,dim,row_eq,col_u);//U
        // v-terms
        // w-terms
        // P-terms
        set_Mat(1./rho,I,ny,B,dim,row_eq,col_P);//1/rho

        // set v-mom equation for A
        row_eq = (4*zi+1)*ny;
        // u-terms
        // v-terms
        set_Mat(-i*m*omega,I,ny,A,dim,row_eq,col_v);//-imw
        set_Mat(i*alpha,U,ny,A,dim,row_eq,col_v);//ialpha*U
        set_Mat((1./Re) * (pow(alpha,2)),I,ny,A,dim,row_eq,col_v);//-Re^-1 - alpha^2
        set_Mat(-1./Re,Dyy,ny,A,dim,row_eq,col_v);//-Re^-1 Dyy
        set_Mat(-1./Re,Dzz,nz,ny,zi,A,dim,row_eq,col_v-4*zi*ny); // -Re^-1 Dzz
        // w-terms
        // P-terms
        set_Mat(1./rho,Dy,ny,A,dim,row_eq,col_P);//1/rho Dy
        // set v-mom equation for B
        // u-terms
        // v-terms
        set_Mat(1.,U,ny,B,dim,row_eq,col_v);//U
        // w-terms
        // P-terms

        // set w-mom equation for A
        row_eq = (4*zi+2)*ny;
        // u-terms
        // v-terms
        // w-terms
        set_Mat(-i*m*omega,I,ny,A,dim,row_eq,col_w);//-imw
        set_Mat(i*alpha,U,ny,A,dim,row_eq,col_w);//ialpha*U
        set_Mat((1./Re) * (pow(alpha,2)),I,ny,A,dim,row_eq,col_w);//-Re^-1 - alpha^2
        set_Mat(-1./Re,Dyy,ny,A,dim,row_eq,col_w);//-Re^-1 - alpha^2
        set_Mat(-1./Re,Dzz,nz,ny,zi,A,dim,row_eq,col_w-4*zi*ny); // -Re^-1 Dzz
        // P-terms
        //set_Mat(1./rho,Dy,ny,A,dim,row_eq,col_P);//1/rho Dy
        set_Mat(1./rho,Dz,nz,ny,zi,A,dim,row_eq,col_P-4*zi*ny); // rho^-1 Dz
        // set w-mom equation for B
        // u-terms
        // v-terms
        // w-terms
        set_Mat(1.,U,ny,B,dim,row_eq,col_w);//U
        // P-terms
        
        // set continuity equation for A
        row_eq = (4*zi+3)*ny;
        // u-terms
        set_Mat(i*alpha,I,ny,A,dim,row_eq,col_u);//ialpha*u
        // v-terms
        set_Mat(1.,Dy,ny,A,dim,row_eq,col_v);//Dy*v
        // w-terms
        set_Mat(1.,Dz,nz,ny,zi,A,dim,row_eq,col_w-4*zi*ny); // Dz*w
        // P-terms
        // set continuity equation for B
        // u-terms
        set_Mat(1.,I,ny,B,dim,row_eq,col_u);//I*u
        // v-terms
        // w-terms
        // P-terms
        
        return 0;
    }
}
