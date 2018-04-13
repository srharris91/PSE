#include "print.hpp"
namespace PSE
{
    PetscErrorCode printScalar( PetscScalar x[], int n,char const name[]){
        PetscErrorCode ierr;
        for(int i=0; i<n; i++) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"  %s[%D] = %g + %g i\n",name,i,(double)PetscRealPart(x[i]),(double)PetscImaginaryPart(x[i]));CHKERRQ(ierr);}
        return ierr;
    }
    PetscErrorCode printVec( Vec &x, int n,char const name[]){
        PetscErrorCode ierr;
        PetscScalar *xa;
        ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
        for(int i=0; i<n; i++) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"  %s[%D] = %g + %g i\n",name,i,(double)PetscRealPart(xa[i]),(double)PetscImaginaryPart(xa[i]));CHKERRQ(ierr);}
        ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);
        return ierr;
    }
}
