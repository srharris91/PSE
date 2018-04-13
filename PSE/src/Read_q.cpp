#include "Read_q.hpp"
#include "print.hpp"

namespace PSE
{
    int Read_q(
            PetscScalar *output,
            int n,
            char const buff[]
            ){
        PetscErrorCode ierr;
        //char buff[]="tofile.dat";
        FILE *latfile;
        latfile=fopen(buff,"r");
        fread(output,sizeof(double),n,latfile);
        fclose(latfile);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nReading in vector\n");CHKERRQ(ierr);
        for (int i=0; i<n; i++){
            ierr = PetscPrintf(PETSC_COMM_WORLD,"x[%i] = %g + %g i\n",
                    i,
                    (double)PetscRealPart(output[i]),
                    (double)PetscImaginaryPart(output[i])); CHKERRQ(ierr);
        }
        //printScalar(output,n);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"done Reading in output\n\n\n");CHKERRQ(ierr);
        return 0;

    }
}
