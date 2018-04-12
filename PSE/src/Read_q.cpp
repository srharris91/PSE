#include "Read_q.hpp"

namespace PSE
{
    int Read_q(
            PetscScalar *RHS_True,
            int n){
        PetscErrorCode ierr;
        char buff[]="tofile.dat";
        FILE *latfile;
        latfile=fopen(buff,"r");
        fread(RHS_True,sizeof(double),n,latfile);
        fclose(latfile);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nReading in RHS_True\n");CHKERRQ(ierr);
        for(int i=0; i<n; i++){
            ierr = PetscPrintf(PETSC_COMM_WORLD,"  RHS_True[%D] = %g \n",i,(double)RHS_True[i]);CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD,"done Reading in RHS_True\n\n\n");CHKERRQ(ierr);
        return 0;

    }
}
