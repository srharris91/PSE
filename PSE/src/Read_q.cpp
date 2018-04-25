#include "Read_q.hpp"
#include "print.hpp"
#include "set_Vec.hpp"
#include "Init_Vec.hpp"

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
        if(rank==0){
            FILE *latfile;
            latfile=fopen(buff,"r");
            fread(read_scalar,sizeof(double),2*n,latfile);
            fclose(latfile);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nReading in vector\n");CHKERRQ(ierr);
            /*
            for (int i=0; i<n; i++){// 2*n because of complex variables
                ierr = PetscPrintf(PETSC_COMM_WORLD,"  x[%i] = %g + %g i\n",
                        i,
                        (double)PetscRealPart(read_scalar[i]),
                        (double)PetscImaginaryPart(read_scalar[i])); CHKERRQ(ierr);
            }
            */
            //printScalar(read_scalar,n);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"done Reading in read_scalar\n\n\n");CHKERRQ(ierr);
            set_Vec(read_scalar,n,output,PETSC_FALSE); // set vector on rank = 0
        }
        set_Vec(output); // assemble matrix on all processors
        return 0;

    }
    /* still need to add this feature of reading in with petsc viewer from binary file 
     PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from vector.dat ...\n");
     PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector.dat",FILE_MODE_READ,&viewer);
     VecCreate(PETSC_COMM_WORLD,&u);
     VecLoad(u,viewer);
     PetscViewerDestroy(&viewer);
     */
    /*
    int Read_q(
            Vec &output,
            int n,
            char const buff[]
            ){
        PetscViewer     viewer;
        PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from %s ...\n",buff);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,buff,FILE_MODE_READ,&viewer);
        //VecCreate(PETSC_COMM_WORLD,&output);
        Init_Vec(output,n);
        VecLoad(output,viewer);
        PetscViewerDestroy(&viewer);
        return 0;
    }
    */

}
