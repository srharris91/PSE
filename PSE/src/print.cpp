#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#include "print.hpp"
namespace PSE
{
    PetscErrorCode printScalar( 
            const PetscScalar x[],
            const PetscInt n,
            char const name[],
            const PetscViewer viewer){
        PetscErrorCode ierr;
        if (n==1){
            ierr = PetscPrintf(PETSC_COMM_WORLD,ANSI_COLOR_CYAN "  %s = %e + %e i\n" ANSI_COLOR_RESET,name,(double)PetscRealPart(x[0]),PetscImaginaryPart(x[0]));CHKERRQ(ierr);
        }
        else{
            PetscScalarView(n,x,viewer);
        }
        return ierr;
    }
    PetscErrorCode printVec( 
            const Vec &x,
            const PetscInt n,
            char const name[]){
        PetscErrorCode ierr;
        PetscScalar *xa;
        ierr = VecGetArray(x,&xa);CHKERRQ(ierr);
        for(PetscInt i=0; i<n; i++) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"  %s[%D] = %g + %g i\n",name,i,(double)PetscRealPart(xa[i]),(double)PetscImaginaryPart(xa[i]));CHKERRQ(ierr);
        }

        ierr = VecRestoreArray(x,&xa);CHKERRQ(ierr);
        return ierr;
    }
    PetscErrorCode printInt( 
            const PetscInt x[],
            const PetscInt n,
            char const name[]){
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
    void printVecView( 
            const Vec &x,
            char const name[],
            const PetscViewerFormat format){
        PetscViewer     viewer;
        PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL,NULL, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,&viewer);
        PetscObjectSetName((PetscObject)viewer,name);
        PetscViewerPushFormat(viewer,format);

        VecView(x,viewer);
        PetscViewerDestroy(&viewer);
    }
    void printMatView( 
            const Mat &A,
            char const name[],
            const PetscViewerFormat format
            ){
        PetscViewer     viewer;
        PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL,NULL, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,&viewer);
        PetscObjectSetName((PetscObject)viewer,name);
        PetscViewerPushFormat(viewer,format);
        MatView(A,viewer);
        PetscViewerDestroy(&viewer);
    }
    void printMatASCII( 
            const Mat &A,
            char const name[],
            const PetscViewerFormat format
            ){
        PetscViewer     viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD,name,&viewer);
        PetscViewerPushFormat(viewer,format);
        MatView(A,viewer);
        PetscViewerDestroy(&viewer);
    }
    void printVecASCII( 
            const Vec &b,
            char const name[],
            const PetscViewerFormat format
            ){
        PetscViewer     viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD,name,&viewer);
        PetscViewerPushFormat(viewer,format);
        VecView(b,viewer);
        PetscViewerDestroy(&viewer);
    }
}
