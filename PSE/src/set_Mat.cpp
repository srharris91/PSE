#include "set_Mat.hpp"
#include "Init_Mat.hpp"
#include "print.hpp"
#include <iostream>

namespace PSE{
    PetscInt set_Mat(
            const PetscScalar* const* Ain,
            const PetscInt &n,
            Mat &A
            ){
        PetscInt Istart, Iend, Ii, col[n],i;
        for(i=0;i<n;i++) col[i]=i;
        PetscErrorCode ierr;

        // initialize A matrix of right size
        Init_Mat(A,n);

        // Currently, all PETSc parallel matrix formats are partitioned by
        // contiguous chunks of rows across the processors.  Determine which
        // rows of the matrix are locally owned.
        ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

        // Set matrix elements (in parallel)

        for (Ii=Istart; Ii<Iend; Ii++){ // in parallel all of the values
            //std::cout<<"setting i="<<Ii<<" n="<<n<<" col[:] = ";
            //for(i=0;i<n;i++) std::cout<<col[i]<<", ";
            //std::cout<<std::endl;
            ierr = MatSetValues(A,1,&Ii,n,col,Ain[Ii],ADD_VALUES);CHKERRQ(ierr);
        }
        //std::cout<<"done setting i="<<Ii<<std::endl;


        // Assemble matrix, using the 2-step process:
        // MatAssemblyBegin(), MatAssemblyEnd()
        // Computations can be done while messages are in transition
        // by placing code between these two statements.
        //ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        //ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

        return 0;
    }
    PetscInt set_Mat(
            const PetscScalar &diag,
            const PetscInt &n,
            Mat &A,
            const PetscInt &k,
            const PetscBool &parallel
            ){
        PetscErrorCode ierr;

        // initialize A matrix of right size
        //Init_Mat(A,n);


        // Set matrix elements (in parallel)

        if (parallel){
            PetscInt Istart, Iend, Ii;
            // Currently, all PETSc parallel matrix formats are partitioned by
            // contiguous chunks of rows across the processors.  Determine which
            // rows of the matrix are locally owned.
            ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
            for (Ii=Istart; Ii<Iend; Ii++){ // in parallel all of the values
                if (Ii+k >= 0 && Ii+k <= n-1){
                    ierr = MatSetValue(A,Ii,(Ii + k),diag,ADD_VALUES);CHKERRQ(ierr);
                }
            }
        }
        else{
            for (PetscInt Ii=0; Ii<n; Ii++){ // in parallel all of the values
                if (Ii+k >= 0 && Ii+k <= n-1){
                    ierr = MatSetValue(A,Ii,(Ii + k),diag,ADD_VALUES);CHKERRQ(ierr);
                }
            }
        }


        return 0;
    }
    PetscInt set_Mat(
            Mat &A  //< Mat to assemble
            ){
        // Assemble matrix, using the 2-step process:
        // MatAssemblyBegin(), MatAssemblyEnd()
        // Computations can be done while messages are in transition
        // by placing code between these two statements.

        PetscErrorCode ierr;
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        return 0;
    }
    PetscInt set_Mat(
            const PetscScalar &value,
            const PetscInt &row,
            const PetscInt &col,
            Mat &A,
            const PetscBool &parallel
            ){
        PetscErrorCode ierr;
        if (parallel){
            PetscInt Istart, Iend, Ii;
            ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
            for (Ii=Istart; Ii<Iend; Ii++){ // in parallel all of the values
                if (Ii == row){
                    ierr = MatSetValue(A,row,col,value,ADD_VALUES);CHKERRQ(ierr);
                }
            }
        }
        else{
            PetscMPIInt rank;
            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            if(rank==0){
                ierr = MatSetValue(A,row,col,value,ADD_VALUES);CHKERRQ(ierr);
            }
        }
        return 0;
    }
    PetscInt set_Mat(
            const PetscScalar &a,
            const Mat &Asub,
            const PetscInt &nsub,
            Mat &A,
            const PetscInt &n,
            const PetscInt &rowoffset,
            const PetscInt &coloffset,
            const InsertMode &addv
            ){
        PetscInt Isubstart,Isubend;
        PetscInt ncols;
        const PetscInt *cols;
        const PetscScalar *vals;

        // get ranges for each matrix
        MatGetOwnershipRange(Asub,&Isubstart,&Isubend);
        //MatGetOwnershipRange(A,&Istart,&Iend);

        //MatGetRow(Asub,0,&ncols,&cols,&vals);
        //MatRestoreRow(Asub,0,&ncols,&cols,&vals);
        for (PetscInt i=Isubstart; i<Isubend; i++){
            MatGetRow(Asub,i,&ncols,&cols,&vals);
            PetscInt offcols[ncols];
            PetscScalar avals[ncols];
            for (PetscInt j=0; j<ncols; j++) {
                offcols[j] = cols[j]+coloffset;
                avals[j] = a*vals[j];
            }
            //std::cout<<"vals [i] = "<<vals[i-Isubstart]<<std::endl;
            //printScalar(vals,ncols);
            set_Mat(avals,i+rowoffset,ncols,offcols,A,addv);
            
            MatRestoreRow(Asub,i,&ncols,&cols,&vals);
        }
        return 0;
    }

    PetscInt set_Mat(
            const PetscScalar* Ain,
            const PetscInt &row,
            const PetscInt &ncols,
            const PetscInt cols[],
            Mat &A,
            const InsertMode &addv
            ){
        PetscErrorCode ierr;
        ierr = MatSetValues(A,1,&row,ncols,cols,Ain,addv);CHKERRQ(ierr);
        return 0;
    }

    PetscInt set_Mat(
            const PetscScalar &a,
            const Mat &Dz,           
            const PetscInt &nz,     
            const PetscInt &ny,    
            const PetscInt &zi,
            Mat &A,               
            const PetscInt &n,   
            const PetscInt &rowoffset,
            const PetscInt &coloffset,
            const InsertMode &addv
            ){
        PetscInt Isubstart,Isubend;
        PetscInt ncols;
        const PetscInt *cols;
        const PetscScalar *vals;

        // get ranges for each matrix
        MatGetOwnershipRange(Dz,&Isubstart,&Isubend);

        for (PetscInt i=Isubstart; i<Isubend; i++){
            if(i==zi){ // if this processor contains the correct zi plane
                MatGetRow(Dz,i,&ncols,&cols,&vals);
                PetscInt offcols[ncols];
                PetscScalar avals[ncols];
                for (PetscInt k=0; k<ny; k++){ // for each row in A
                    for (PetscInt j=0; j<ncols; j++) {
                        offcols[j] = 4*ny*cols[j]+coloffset+k;
                        avals[j] = a*vals[j];
                    }
                    //std::cout<<"vals [i] = "<<vals[i-Isubstart]<<std::endl;
                    //printScalar(vals,ncols);
                    set_Mat(avals,i+rowoffset+k,ncols,offcols,A,addv);
                }

                MatRestoreRow(Dz,i,&ncols,&cols,&vals);
            }
        }
        return 0;
    }
}

