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
        Init_Mat(output,n);
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
        PetscScalar smax = s[Nm1-1];

        // solve for Coefficients given stencil
        Vec CoeffsMain;
        //Init_Vec(CoeffsMain,Nm1);
        get_D_Coeffs(s,Nm1,CoeffsMain,d);
        //printVecView(Coeffs,Nm1);

        //set diagonals
        Vec sVec;
        Init_Vec(sVec,Nm1);
        set_Vec(s,Nm1,sVec);
        //printVecView(sVec,Nm1);
        // set D matrix
        VecGetLocalSize(CoeffsMain,&nlocal);
        VecGetArray(CoeffsMain,&array);
        VecGetArray(sVec,&sarray);
        // set main diagonals
        PetscBool parallel=PETSC_FALSE;
        for (i=0; i<nlocal; i++) { set_Mat(array[i]/denom,n,output,(PetscInt)PetscRealPart(sarray[i]),parallel); }
        // set periodic diagonals
        if (periodic){
            for (i=0; i<nlocal; i++) { 
                if ((PetscInt)PetscRealPart(sarray[i]) > 0){
                    set_Mat(array[i]/denom,n,output,(PetscInt)PetscRealPart(sarray[i])-n,parallel);
                }
                else if ((PetscInt)PetscRealPart(sarray[i]) < 0){
                    set_Mat(array[i]/denom,n,output,(PetscInt)PetscRealPart(sarray[i])+n,parallel);
                }
            }
        }
        set_Mat(output); // do nothing, but assemble
        VecRestoreArray(CoeffsMain,&array);
        VecRestoreArray(sVec,&sarray);

        // set BC so we don't go out of range on top and bottom
        if (!periodic){
            //if (d % 2 !=0) Nm1-=1;// need one less pt in shifted diff if odd derivative (vs central)
            Vec Coeffsbottom[(PetscInt)PetscRealPart(smax)];
            Vec Coeffstop[(PetscInt)PetscRealPart(smax)];
            for (i=0; i<PetscRealPart(smax); i++){
                // bottom BC
                if (reduce_wall_order){
                    PetscScalar *s = new PetscScalar[Nm1]; // resize
                    for (j=0; j<Nm1; j++) s[j] = j - i; // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,Nm1,Coeffsbottom[i],d);
                    // set values in matrix
                    set_Vec(s,Nm1,sVec);
                    //printVecView(sVec,Nm1);
                }
                else{
                    PetscScalar *s = new PetscScalar[N]; // resize
                    for (j=0; j<N; j++) s[j] = j - i; // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,N,Coeffsbottom[i],d);
                    // set values in matrix
                    set_Vec(s,N,sVec);
                }
                // set D matrix
                VecGetLocalSize(Coeffsbottom[i],&nlocal);
                VecGetArray(Coeffsbottom[i],&array);
                VecGetArray(sVec,&sarray);
                MatZeroRows(output,1,&i,0.,0,0);
                for (j=0; j<nlocal; j++) { 
                    MatSetValue(output,i,(PetscInt)PetscRealPart(sarray[j])+i,array[j]/denom,INSERT_VALUES);
                }
                set_Mat(output); // do nothing, but assemble
                VecRestoreArray(Coeffsbottom[i],&array);
                ierr = VecDestroy(&Coeffsbottom[i]);CHKERRQ(ierr);
                VecRestoreArray(sVec,&sarray);
 
                // top BC
                if (reduce_wall_order){
                    PetscScalar *s = new PetscScalar[Nm1]; // resize
                    for (j=0; j<Nm1; j++) s[j] = -(j - i); // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,Nm1,Coeffstop[i],d);
                    // set values in matrix
                    set_Vec(s,Nm1,sVec);
                    //printVecView(sVec,Nm1);
                }
                else{
                    PetscScalar *s = new PetscScalar[N]; // resize
                    for (j=0; j<N; j++) s[j] = -(j - i); // stencil over central diff of order
                    // solve for Coefficients given stencil
                    get_D_Coeffs(s,N,Coeffstop[i],d);
                    // set values in matrix
                    set_Vec(s,N,sVec);
                }
                // set D matrix
                VecGetLocalSize(Coeffstop[i],&nlocal);
                VecGetArray(Coeffstop[i],&array);
                VecGetArray(sVec,&sarray);
                PetscInt rowi = n-i-1;
                MatZeroRows(output,1,&rowi,0.,0,0);
                for (j=0; j<nlocal; j++) { 
                    MatSetValue(output,rowi,rowi+(PetscInt)PetscRealPart(sarray[j]),array[j]/denom,INSERT_VALUES);
                }
                set_Mat(output); // do nothing, but assemble
                VecRestoreArray(Coeffstop[i],&array);
                VecRestoreArray(sVec,&sarray);
                ierr = VecDestroy(&Coeffstop[i]);CHKERRQ(ierr);
            }
        }
        //printVecView(Coeffs,Nm1);
        //MatDuplicate(output,MAT_COPY_VALUES,&output);
        

        delete[] s;
        ierr = VecDestroy(&CoeffsMain);CHKERRQ(ierr);
        ierr = VecDestroy(&sVec);CHKERRQ(ierr);
        //ierr = MatDestroy(&D); CHKERRQ(ierr);

        return 0;
    }
}
