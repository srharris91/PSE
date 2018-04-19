#include "print.hpp"
#include <iostream>
namespace PSE
{
    PetscErrorCode printScalar( const PetscScalar x[], const PetscInt n,char const name[],PetscViewer viewer){
        PetscErrorCode ierr;
        if (n==1){
            ierr = PetscPrintf(PETSC_COMM_WORLD,"  %s = %g \n",name,(double)PetscRealPart(x[0]));CHKERRQ(ierr);
        }
        else{
            PetscScalarView(n,x,viewer);
        }
        return ierr;
    }
    PetscErrorCode printVec( Vec &x, PetscInt n,char const name[]){
        PetscErrorCode ierr;
        PetscScalar *xa;
        ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
        for(PetscInt i=0; i<n; i++) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"  %s[%D] = %g + %g i\n",name,i,(double)PetscRealPart(xa[i]),(double)PetscImaginaryPart(xa[i]));CHKERRQ(ierr);
            //std::cout<<"     xa["<<i<<"]= "<<xa[i]<<std::endl;
        }

        ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);
        return ierr;
    }
    PetscErrorCode printInt( PetscInt x[], PetscInt n,char const name[]){
        PetscErrorCode ierr;
        if (n==1){
            ierr = PetscPrintf(PETSC_COMM_WORLD,"  %s = %i \n",name,(PetscInt)x[0]);CHKERRQ(ierr);
        }
        else{
            PetscErrorCode ierr;
            for(PetscInt i=0; i<n; i++) {
                ierr = PetscPrintf(PETSC_COMM_WORLD,"  %s[%D] = %i \n",name,i,(PetscInt)x[i]);CHKERRQ(ierr);
            }
        }
        return ierr;
    }
    void printVecView( Vec &x, PetscInt n,char const name[],PetscViewerFormat format){
        PetscViewer     viewer;
        PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL,NULL, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,&viewer);
        PetscObjectSetName((PetscObject)viewer,name);
        PetscViewerPushFormat(viewer,format);

        VecView(x,viewer);
    }
    void printMatView( 
            const Mat &A,
            const PetscInt n,
            char const name[],
            const PetscViewerFormat format
            ){
        PetscViewer     viewer;
        PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL,NULL, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,&viewer);
        PetscObjectSetName((PetscObject)viewer,name);
        PetscViewerPushFormat(viewer,format);

        MatView(A,viewer);
    }
}
