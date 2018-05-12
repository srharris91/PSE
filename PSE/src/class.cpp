#include "class.hpp"
#include "Init_Mat.hpp"
#include "Init_Vec.hpp"
#include "set_Mat.hpp"

namespace PSE
{
    PSE::PSE(){
        // set default parameters
        Re = 2000.;
        rho = 1.;
        m=1.;
        omega=0.3;
        ny=101;
        nz=6;
        hx=2.5;
        xsteps=15;
        dim=ny*nz*4;
        order=4;

        // initialize matrices
        init();
    }
    void PSE::init(){
        if (flag_init==PETSC_FALSE){
            // initialize Matrices in this class
            // initialize and set I
            if (flag_set_I==PETSC_FALSE){
                Init_Mat(I,ny);
                set_Mat(1.,ny,I); // set diag of 1
                set_Mat(I); // assemble
                flag_set_I=PETSC_TRUE;
            }
            Init_Mat(A,dim);
            Init_Mat(B,dim);

            // initialize Vectors in this class
            Init_Vec(q,dim);
            Init_Vec(qp1,dim);
            Init_Vec(b,dim);
            
            // initialize dynamic arrays
            y=new PetscScalar[ny];
            z=new PetscScalar[nz];
            
            
            flag_init=PETSC_TRUE;
        }


    }
    PSE::~PSE(){
        // remove PetscScalar dynamic arrays
        delete[] y;
        delete[] z;
    }

    PetscInt PSE::destroy(){
        ierr = MatDestroy(&Dy); CHKERRQ(ierr);
        ierr = MatDestroy(&Dyy);CHKERRQ(ierr);
        ierr = MatDestroy(&Dz); CHKERRQ(ierr);
        ierr = MatDestroy(&Dzz);CHKERRQ(ierr);
        ierr = MatDestroy(&I);  CHKERRQ(ierr);
        ierr = MatDestroy(&U);  CHKERRQ(ierr);
        ierr = MatDestroy(&Uy); CHKERRQ(ierr);
        ierr = MatDestroy(&A);  CHKERRQ(ierr);
        ierr = MatDestroy(&B);  CHKERRQ(ierr);
        ierr = VecDestroy(&q);  CHKERRQ(ierr);
        ierr = VecDestroy(&qp1);CHKERRQ(ierr);
        ierr = VecDestroy(&b)  ;CHKERRQ(ierr);
    }
}
