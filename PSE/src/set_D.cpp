#include <stdexcept>
#include <string>
#include "set_D.hpp"
#include "get_D_Coeffs.hpp"
#include "Init_Vec.hpp"
#include "Init_Mat.hpp"
#include "print.hpp"
#include <iostream>

namespace PSE
{
    PetscInt set_D(
            const PetscScalar y[],
            const PetscInt &n,   
            Mat &output,         
            const PetscInt &order,
            const PetscInt &d ,   
            const PetscBool &reduce_wall_order, 
            const PetscBool &output_full, 
            const PetscBool &periodic
            ){
        PetscErrorCode ierr;
        double h=y[1]-y[0]; // assume uniform spacing
        PetscInt N=order+d; // number of pts needed for order of accuracy
        if (N<n) {
            std::runtime_error(
                    "You need more points in your domain, you need "
                    + std::to_string(N) 
                    + "pts and you only gave " 
                    + std::to_string(n));
        }
        PetscInt Nm1 = N-1; // how many pts needed if using central diff
        if (d % 2 !=0) Nm1+=1;// need one more pt in central diff if odd derivative
        PetscScalar *s = new PetscScalar[Nm1];
        for (int i=0; i<Nm1; i++) s[i] = i - (int((Nm1-1)/2)); // stencil over central diff of order
        for (int i=0; i<Nm1; i++) std::cout<<"s["<<i<<"] = "<<s[i]<<std::endl;; // stencil over central diff of order
        PetscScalar smax = s[Nm1];

        Vec Coeffs;
        Init_Vec(Coeffs,Nm1);
        get_D_Coeffs(s,Nm1,Coeffs);
        //printVecView(Coeffs,Nm1);


        /*
        // create A dynamic 2d array of correct square shape
        PetscScalar **A=new PetscScalar*[n];
        for (int i=0; i<n; i++) A[i]=new PetscScalar[n];
        // create b dynamic 1d matrix
        PetscScalar *b=new PetscScalar[n];

        // set A matrix
        for (int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                A[i][j] = std::pow(s[j],i);
            }
        }
        // set b vector
        for (int i=0; i<n; i++) b[i]=0.;
        b[d]=factorial(d);

        // solve Ax=b problem in parallel using PSE::Ax_b
        ierr = Ax_b(A,output,b,n); CHKERRQ(ierr);
        //printVecView(output,n);
        



        // delete A and b matrix
        for (int i=0; i<n; i++) delete[] A[i];
        delete[] A;
        delete[] b;
        */
        delete[] s;
        ierr = VecDestroy(&Coeffs);CHKERRQ(ierr);

        return 0;
    }
}
