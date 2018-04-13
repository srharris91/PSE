#include "Read_q.hpp"
#include "print.hpp"
#include <mpi.h>
#include <iostream>

namespace PSE
{
    int Read_q(
            Vec &output,
            int n,
            char const buff[]
            ){
        PetscErrorCode ierr;
        PetscScalar read_scalar[n];
        for(int i=0; i<n; i++) read_scalar[i]=0;
        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        //char buff[]="tofile.dat";
        if(rank==0){
            FILE *latfile;
            latfile=fopen(buff,"r");
            fread(read_scalar,sizeof(double),n,latfile);
            fclose(latfile);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nReading in vector\n");CHKERRQ(ierr);
            for (int i=0; i<n; i++){
                ierr = PetscPrintf(PETSC_COMM_WORLD,"x[%i] = %g + %g i\n",
                        i,
                        (double)PetscRealPart(read_scalar[i]),
                        (double)PetscImaginaryPart(read_scalar[i])); CHKERRQ(ierr);
            }
            //printScalar(read_scalar,n);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"done Reading in read_scalar\n\n\n");CHKERRQ(ierr);
        }
        for (int i=0; i<n; i++){
            ierr = VecSetValues(output,1,&i,&read_scalar[i],ADD_VALUES);CHKERRQ(ierr);
        }
        ierr = VecAssemblyBegin(output); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(output); CHKERRQ(ierr);
        return 0;

    }
}
