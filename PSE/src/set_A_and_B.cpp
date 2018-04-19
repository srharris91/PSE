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
#include "base_flow.hpp"
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
            const PetscScalar &Re,
            const PetscScalar &rho,
            const PetscScalar &alpha,
            const PetscScalar &m,
            const PetscScalar &omega,
            const PetscInt &order,
            const PetscBool &reduce_wall_order
            ){
        PetscErrorCode ierr;
        PetscInt dim=ny*nz*4; // dimension of A and B square matrices
        PetscScalar i = 1.*PETSC_i;
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
        set_D(z,nz,Dz,order,1);
        set_D(z,nz,Dzz,order,2);
        set_D(y,ny,Dyy,order,2);
        
        Init_Mat(A,dim);
        Init_Mat(B,dim);

        PetscInt row_eq;
        PetscInt col_u=0,
                 col_v=ny,
                 col_w=2*ny,
                 col_P=3*ny;
        // set u-mom equation for A
        row_eq=0;
        // u-terms
        set_Mat(-i*m*omega,I,ny,A,dim);//-imw
        set_Mat(i*alpha,U,ny,A,dim);//ialpha*U
        //set_Mat(1.,Uy,ny,A,dim);//Ux
        set_Mat((-1./Re) * (-pow(alpha,2)),I,ny,A,dim);//Re^-1 * alpha^2
        set_Mat(-1./Re,Dyy,ny,A,dim);// -Re^-1 Dyy
        set_Mat(-1./Re,Dzz,nz,ny,0,A,dim); // -Re^-1 Dzz
        // v-terms
        set_Mat(1.,Uy,ny,A,dim,row_eq,col_v);//Uy
        //w-terms
        //P-terms
        set_Mat((1./rho) * i*alpha,I,ny,A,dim,row_eq,col_P);//imw
        // set u-mom equation for B
        set_Mat(1.,U,ny,B,dim);//ialpha*U
        set_Mat(1./rho,I,ny,B,dim,row_eq,col_P);//ialpha*U

        // set v-mom equation for A
        row_eq = ny;
        // u-terms
        // v-terms
        set_Mat(-i*m*omega,I,ny,A,dim,row_eq,col_v);//-imw
        set_Mat(i*alpha,U,ny,A,dim,row_eq,col_v);//ialpha*U
        set_Mat((-1./Re) * (-pow(alpha,2)),I,ny,A,dim,row_eq,col_v);//-Re^-1 - alpha^2
        set_Mat(-1./Re,Dyy,ny,A,dim,row_eq,col_v);//-Re^-1 - alpha^2
        set_Mat(-1./Re,Dzz,nz,ny,0,A,dim,row_eq,col_v); // -Re^-1 Dzz
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
        row_eq = 2*ny;
        // u-terms
        // v-terms
        // w-terms
        set_Mat(-i*m*omega,I,ny,B,dim,row_eq,col_w);//-imw
        set_Mat(i*alpha,U,ny,A,dim,row_eq,col_w);//ialpha*U
        set_Mat((-1./Re) * (-pow(alpha,2)),I,ny,A,dim,row_eq,col_w);//-Re^-1 - alpha^2
        set_Mat(-1./Re,Dyy,ny,A,dim,row_eq,col_w);//-Re^-1 - alpha^2
        set_Mat(-1./Re,Dzz,nz,ny,0,A,dim,row_eq,col_w); // -Re^-1 Dzz
        // P-terms
        //set_Mat(1./rho,Dy,ny,A,dim,row_eq,col_P);//1/rho Dy
        set_Mat(1./rho,Dz,nz,ny,0,A,dim,row_eq,col_P); // rho^-1 Dz
        // set w-mom equation for B
        // u-terms
        // v-terms
        // w-terms
        set_Mat(1.,U,ny,B,dim,row_eq,col_w);//U
        // P-terms
        
        // set continuity equation for A
        row_eq = 3*ny;
        // u-terms
        set_Mat(-i*alpha,I,ny,A,dim,row_eq,col_u);//ialpha*u
        // v-terms
        set_Mat(1.,Dy,ny,A,dim,row_eq,col_v);//Dy*v
        // w-terms
        set_Mat(1.,Dz,nz,ny,0,A,dim,row_eq,col_w); // Dz*w
        // P-terms
        // set continuity equation for B
        // u-terms
        set_Mat(1.,I,ny,B,dim,row_eq,col_u);//I*u
        // v-terms
        // w-terms
        // P-terms
        

        //set_Mat(1.,I,ny,A,dim);//imw
        
       
      
        //set_Mat(1.,Dy,ny,B,dim);
        //set_Mat(1.,Dz,nz,ny,0,B,dim);
        //set_Mat(1.,Dzz,nz,ny,0,B,dim);

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
