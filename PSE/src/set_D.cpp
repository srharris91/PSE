#include <stdexcept>
#include <string>
#include <cmath>
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
    PetscInt set_D(
            const PetscScalar y[],
            const PetscInt &n,   
            Mat &output,         
            const PetscInt &order,
            const PetscInt &d ,   
            const PetscBool &periodic,
            const PetscBool &reduce_wall_order
            ){
        //Mat D;
        //Init_Mat(D,n);
        PetscErrorCode ierr;
        PetscScalar *array,*sarray;
        PetscInt i, j, nlocal;
        PetscScalar h=y[1]-y[0]; // assume uniform spacing
        PetscScalar denom=pow(h,d); // denominator h value
        PetscInt N=order+d; // number of pts needed for order of accuracy
        if (N>n) {
            throw std::runtime_error(
                    "You need more points in your domain, you need "
                    + std::to_string(N) 
                    + "pts and you only gave " 
                    + std::to_string(n));
        }
        PetscInt Nm1 = N-1; // how many pts needed if using central diff
        if (d % 2 !=0) Nm1+=1;// need one more pt in central diff if odd derivative
        PetscScalar *s = new PetscScalar[Nm1];
        for (i=0; i<Nm1; i++) s[i] = i - (int((Nm1-1)/2)); // stencil over central diff of order
        //for (int i=0; i<Nm1; i++) std::cout<<"s["<<i<<"] = "<<s[i]<<std::endl;; // stencil over central diff of order
        PetscScalar smax = s[Nm1-1];

        // solve for Coefficients given stencil
        Vec Coeffs;
        Init_Vec(Coeffs,Nm1);
        get_D_Coeffs(s,Nm1,Coeffs);
        //printVecView(Coeffs,Nm1);

        //set diagonals
        Vec sVec;
        Init_Vec(sVec,Nm1);
        set_Vec(s,Nm1,sVec);
        //printVecView(sVec,Nm1);
        // set D matrix
        VecGetLocalSize(Coeffs,&nlocal);
        VecGetArray(Coeffs,&array);
        VecGetArray(sVec,&sarray);
        // set main diagonals
        PetscBool parallel=PETSC_FALSE;
        for (i=0; i<nlocal; i++) { set_Mat(array[i]/denom,n,output,sarray[i],parallel); }
        // set periodic diagonals
        if (periodic){
            for (i=0; i<nlocal; i++) { 
                if (sarray[i] > 0){
                    set_Mat(array[i]/denom,n,output,sarray[i]-n,parallel);
                }
                else if (sarray[i] < 0){
                    set_Mat(array[i]/denom,n,output,sarray[i]+n,parallel);
                }
            }
        }
        set_Mat(output); // do nothing, but assemble
        VecRestoreArray(Coeffs,&array);
        VecRestoreArray(sVec,&array);

        // set BC so we don't go out of range on top and bottom
        if (!periodic){
            if (d % 2 !=0) Nm1-=1;// need one less pt in central diff if odd derivative (vs central)
            for (i=0; i<smax; i++){
                // bottom BC
                if (reduce_wall_order){
                    PetscScalar *s = new PetscScalar[Nm1]; // resize
                    for (j=0; j<Nm1; j++) s[j] = j - i; // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,Nm1,Coeffs);
                    // set values in matrix
                    set_Vec(s,Nm1,sVec);
                    //printVecView(sVec,Nm1);
                }
                else{
                    PetscScalar *s = new PetscScalar[N]; // resize
                    for (j=0; j<N; j++) s[j] = j - i; // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,N,Coeffs);
                    // set values in matrix
                    set_Vec(s,N,sVec);
                }
                // set D matrix
                VecGetLocalSize(Coeffs,&nlocal);
                VecGetArray(Coeffs,&array);
                VecGetArray(sVec,&sarray);
                for (j=0; j<nlocal; j++) { 
                    MatSetValue(output,i,sarray[j]+i,array[j]/denom,INSERT_VALUES);
                    //std::cout<<" row = "<<i<<" col = "<<sarray[j]+i<<" value = "<<array[j]/denom<<std::endl;
                }
                set_Mat(output); // do nothing, but assemble
                VecRestoreArray(Coeffs,&array);
                VecRestoreArray(sVec,&array);
 
                // top BC
                if (reduce_wall_order){
                    PetscScalar *s = new PetscScalar[Nm1]; // resize
                    for (j=0; j<Nm1; j++) s[j] = -(j - i); // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,Nm1,Coeffs);
                    // set values in matrix
                    set_Vec(s,Nm1,sVec);
                    //printVecView(sVec,Nm1);
                }
                else{
                    PetscScalar *s = new PetscScalar[N]; // resize
                    for (j=0; j<N; j++) s[j] = -(j - i); // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,N,Coeffs);
                    // set values in matrix
                    set_Vec(s,N,sVec);
                }
                // set D matrix
                VecGetLocalSize(Coeffs,&nlocal);
                VecGetArray(Coeffs,&array);
                VecGetArray(sVec,&sarray);
                for (j=0; j<nlocal; j++) { 
                    MatSetValue(output,n-i-1,n+sarray[j]-i-1,array[j]/denom,INSERT_VALUES);
                    //std::cout<<" row = "<<n-i-1<<" col = "<<n+sarray[j]-i-1<<" value = "<<array[j]/denom<<std::endl;
                }
                set_Mat(output); // do nothing, but assemble
                VecRestoreArray(Coeffs,&array);
                VecRestoreArray(sVec,&array);
            }
        }
        //printVecView(Coeffs,Nm1);
        //MatDuplicate(output,MAT_COPY_VALUES,&output);
        

        delete[] s;
        ierr = VecDestroy(&Coeffs);CHKERRQ(ierr);
        ierr = VecDestroy(&sVec);CHKERRQ(ierr);
        //ierr = MatDestroy(&D); CHKERRQ(ierr);

        return 0;
    }
}
