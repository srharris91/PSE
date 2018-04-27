#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#include "Read_q.hpp"
#include "print.hpp"
#include "set_Vec.hpp"
#include "Init_Vec.hpp"
#include <string>

namespace PSE
{
    int Read_q(
            Vec &evec,
            PetscScalar y[],
            const int ny,
            PetscScalar z[],
            const int nz,
            PetscScalar &alpha,
            char const buff[]
            ){
        std::string fevec(buff), fy(buff), fz(buff), falpha(buff);
        Read_q(evec,ny*nz*4,fevec.append("_evec.dat").c_str());
        Read_q(y,ny,fy.append("_y.dat").c_str());
        Read_q(z,nz,fz.append("_z.dat").c_str());
        Read_q(alpha,falpha.append("_eig.dat").c_str());

        return 0;
    }
    int Read_q(
            Vec &output,
            const int n,
            char const buff[]
            ){
        PetscErrorCode ierr;
        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if(rank==0){
            PetscScalar *read_scalar=new PetscScalar[n];
            for(int i=0; i<n; i++) read_scalar[i]=0;
            FILE *latfile;
            latfile=fopen(buff,"r");
            fread(read_scalar,sizeof(double),2*n,latfile);
            fclose(latfile);
            ierr = PetscPrintf(PETSC_COMM_WORLD,ANSI_COLOR_GREEN "\nReading in %s\n" ANSI_COLOR_RESET,buff);CHKERRQ(ierr);
            /*
            for (int i=0; i<n; i++){// 2*n because of complex variables
                ierr = PetscPrintf(PETSC_COMM_WORLD,"  x[%i] = %g + %g i\n",
                        i,
                        (double)PetscRealPart(read_scalar[i]),
                        (double)PetscImaginaryPart(read_scalar[i])); CHKERRQ(ierr);
            }
            */
            //printScalar(read_scalar,n);
            set_Vec(read_scalar,n,output,PETSC_FALSE); // set vector on rank = 0
            delete[] read_scalar;
        }
        set_Vec(output); // assemble matrix on all processors
        return 0;

    }
    int Read_q(
            PetscScalar output[],
            const int n,
            char const buff[]
            ){
        PetscErrorCode ierr;

        FILE *latfile;
        latfile=fopen(buff,"r");
        fread(output,sizeof(double),2*n,latfile); // 2n because of complex numbers
        fclose(latfile);
        ierr = PetscPrintf(PETSC_COMM_WORLD,ANSI_COLOR_GREEN "\nReading in %s\n" ANSI_COLOR_RESET,buff);CHKERRQ(ierr);
        return 0;

    }
    int Read_q(
            PetscScalar &output,
            char const buff[]
            ){
        PetscErrorCode ierr;

        PetscScalar read_scalar[2];
        FILE *latfile;
        latfile=fopen(buff,"r");
        fread(read_scalar,sizeof(double),2,latfile); // 2n because of complex numbers
        fclose(latfile);
        output = read_scalar[0] + read_scalar[1]*PETSC_i;
        //PetscRealPart(output) = read_scalar[0];
        //PetscImaginaryPart(output) = read_scalar[1];
        ierr = PetscPrintf(PETSC_COMM_WORLD,ANSI_COLOR_GREEN "\nReading in %s\n" ANSI_COLOR_RESET,buff);CHKERRQ(ierr);
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
