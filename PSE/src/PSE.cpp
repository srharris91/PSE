//#define ANSI_COLOR_RED     "\x1b[31m"
//#define ANSI_COLOR_GREEN   "\x1b[32m"
//#define ANSI_COLOR_YELLOW  "\x1b[33m"
//#define ANSI_COLOR_BLUE    "\x1b[34m"
//#define ANSI_COLOR_MAGENTA "\x1b[35m"
//#define ANSI_COLOR_CYAN    "\x1b[36m"
//#define ANSI_COLOR_RESET   "\x1b[0m"

static char help[] = "Solves a linear system in parallel with KSP.\n\n";

/** @file
   \concepts KSP^solving the PSE equations
   \concepts complex numbers;
   \concepts Parabolic Stability Equations for Transition modeling

   \description Solves a complex linear system in parallel with KSP.

   The model problem:
      Solve PSE equation in the channel flow
      Dirichlet b.c.'s on all sides
      Use the 3-D, finite difference stencil.

   Compiling the code:
      This code uses the complex numbers version of PETSc, so configure
      must be run to enable this

  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include "PSE.hpp"

int main(int argc,char **args){
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
    // test Ax=b solver twice
    if (0){
        // init some variables
        PetscInt       n = 6;
        Vec    x,x2; 
        PetscScalar    Atest[n][6] = {
            {1,1,-2,1,3,-4},
            {2,-1,1,2,1,-3},
            {1,3,-3,-1,2,1},
            {5,2,-1,-1,2,1},
            {-3,-1,2,3,1,3},
            {4,3,1,-6,-3,-2} };
        PetscScalar **Atest2 = new PetscScalar*[n];
        for (int i=0; i<n; i++) Atest2[i] = new PetscScalar[n];

        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                Atest2[i][j] = Atest[i][j];
            }
        }
        PetscScalar b_nondynamic[] = {7,20,-15,-3,16,-27};
        PetscScalar    *btest = new PetscScalar[n];
        for (int i=0; i<n; i++) btest[i] = b_nondynamic[i];

        // initialize Petsc
        PSE::Init_Vec(x,n);

        // solve Ax=b problem
        ierr = PSE::Ax_b(Atest2,x,btest,n); CHKERRQ(ierr);
        // output solutions
        PSE::printVecView(x);

        // solve Ax=b problem
        PSE::Init_Vec(x2,n);
        ierr = PSE::Ax_b(Atest2,x2,btest,n); CHKERRQ(ierr);
        PSE::printVecView(x2,"x2");

        // finalize
        for (int i=0; i<n; i++) delete[] Atest2[i];
        delete[] Atest2;
        delete[] btest;
        ierr = VecDestroy(&x);CHKERRQ(ierr);
        ierr = VecDestroy(&x2);CHKERRQ(ierr);
    }
    // test Read_q
    if (0){
        PetscInt n=6;
        Vec q;
        PSE::Init_Vec(q,n);
        PSE::Read_q(q,n);
        PSE::printVecView(q,"After reading vector = ");
        //PSE::printVecASCII(q,"q_after_reading_vector.txt");
        //PSE::printVecASCII(q,"afterreadingvector.txt");
        ierr = VecDestroy(&q);CHKERRQ(ierr);
    }
    // test get_D_Coeffs
    if (0){
        Vec x;
        //PetscScalar s[]={-3,-2,-1,0};
        PetscInt n=4;
        PetscScalar *s=new PetscScalar[n];
        s[0] = -3;
        s[1] = -2;
        s[2] = -1;
        s[3] = 0;
        Vec sVec;
        PSE::Init_Vec(sVec,n);
        PSE::set_Vec(s,n,sVec);
        PSE::printVecView(sVec);

        PSE::Init_Vec(x,n);


        PSE::get_D_Coeffs(s,n,x);
        PSE::printVecView(x,"after");

        ierr = VecDestroy(&x);CHKERRQ(ierr);
        delete[] s;
    }
    // set D derivative operators
    if(0){
        Mat Dyyp,Dy;
        PetscInt n=10;
        PetscScalar y[n];
        for (int i=0; i<n; i++) y[i] = i;
        // set periodic
        PSE::set_D(y,n,Dyyp,4,1,PETSC_FALSE);
        PSE::printMatView(Dyyp);
        // set non-periodic
        PSE::set_D(y,n,Dy,4,2);
        PSE::printMatView(Dy);


        ierr = MatDestroy(&Dy);CHKERRQ(ierr);
        ierr = MatDestroy(&Dyyp);CHKERRQ(ierr);
    }
    if(0){ // compare against matmult and Ax_b
        // init
        Mat A,B;
        Vec b,q,qp1,qexact,bexact;
        PetscScalar none=-1.0;
        PetscReal norm;
        PetscInt ny=10,nz=6;
        PetscScalar y[ny],z[nz];
        PetscScalar pfive=0.5;
        PetscScalar hx=0.25;
        PetscInt dim=ny*nz*4;
        for (PetscInt i=-ny/2; i<ny/2; i++) y[i] = ((PetscScalar)i) / ((PetscScalar) ny);
        for (PetscInt i=0; i<nz; i++) z[i] = ((PetscScalar)i) / ((PetscScalar) nz);
        // set A,b with q=0.5
        PSE::Init_Vec(q,dim);
        //PSE::Init_Vec(qp1,dim);
        ierr = VecSet(q,pfive);CHKERRQ(ierr);
        PSE::set_Vec(q); // assemble
        PSE::set_A_and_B(y,ny,z,nz,A,B,2000.,1.,1.,1.,1.);
        PSE::set_Euler_Advance(hx,A,B);
        PSE::Init_Vec(b,dim);
        PSE::set_b(B,q,b);
        // set BCs
        PSE::set_BCs(A,b,ny,nz);
        // set exact q and b
        PSE::Init_Vec(qexact,dim);
        PSE::Init_Vec(bexact,dim);
        ierr = VecSet(qexact,pfive);CHKERRQ(ierr);
        PSE::set_Vec(qexact); // assemble
        MatMult(A,qexact,bexact);
        PSE::set_Vec(bexact);
        // view A,B
        PSE::printMatView(A);
        //PSE::printMatASCII(A,"printMatASCII_dense.txt",PETSC_VIEWER_ASCII_DENSE);
        //PSE::printMatASCII(A);
        // solve Ax=b
        PSE::Ax_b(A,qp1,bexact,dim);
        PSE::printVecView(qp1);
        //PSE::printVecASCII(qp1,"printVecqp1.txt");
        //PSE::printVecASCII(qp1,"printVecqp1_dense.txt",PETSC_VIEWER_ASCII_DENSE);
        PSE::printVecView(qexact);
        // check error
        ierr = VecAXPY(qp1,none,qexact);CHKERRQ(ierr);
        ierr = VecNorm(qp1,NORM_2,&norm);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g \n\n",(double)norm);CHKERRQ(ierr);
        // free memory
        ierr = MatDestroy(&A);CHKERRQ(ierr);
        ierr = MatDestroy(&B);CHKERRQ(ierr);
        ierr = VecDestroy(&b);CHKERRQ(ierr);
        ierr = VecDestroy(&bexact);CHKERRQ(ierr);
        ierr = VecDestroy(&q);CHKERRQ(ierr);
        ierr = VecDestroy(&qp1);CHKERRQ(ierr);
        ierr = VecDestroy(&qexact);CHKERRQ(ierr);
    }
    if(0){ // read in eig-value, eig-vectors, y and z
        // init
        Vec q;
        PetscScalar alpha;
        PetscInt ny=201,nz=6;
        PetscScalar y[ny],z[nz];
        PetscInt dim=ny*nz*4;

        // linspace y,z
        PSE::Init_Vec(q,dim);
        if (0){ // create y and z from linspace
            PetscScalar ay=-1.,by=1.;   // [a,b] for y vector
            PetscScalar az=0.,bz=1.;    // [a,b] for z vector
            PetscScalar dy = (by-ay)/((PetscScalar)ny - 1.);
            PetscScalar dz = (bz-az)/((PetscScalar)nz - 1.);
            PetscScalar val;
            PetscInt i;
            for (i=0,val=ay; i<ny; i++,val+=dy) y[i] = val;
            for (i=0,val=az; i<nz; i++,val+=dz) z[i] = val;
        }
        else{// read from files (including y and z
            PSE::Read_q(q,y,ny,z,nz,alpha);
        }
        //PetscPrintf(PETSC_COMM_WORLD,"alpha = %g + %g i",alpha.real(),alpha.imag());
        PetscPrintf(PETSC_COMM_WORLD,"\noutput:\n");
        PSE::printScalar(&alpha);
        PSE::printScalar(y,ny);
        PSE::printScalar(z,nz);
        PSE::printVecView(q);
        //PSE::printVecASCII(q);
        // free memory
        //ierr = MatDestroy(&A);CHKERRQ(ierr);
        ierr = VecDestroy(&q);CHKERRQ(ierr);
    }
    if(1){ // advance q one step, and check growth
        // init
        Mat A,B;
        Vec b,q,qp1;
        PetscScalar Re=6000.,rho=1.,alpha,m=1.,omega=0.27;
        PetscScalar none=-1.0;
        PetscReal norm;
        PetscInt ny=251,nz=6;
        PetscScalar y[ny],z[nz];
        PetscScalar hx=0.001;
        PetscInt dim=ny*nz*4;
        PSE::Init_Vec(q,dim);
        // read in q,y,z,alpha from binary files
        PetscPrintf(PETSC_COMM_WORLD,"\nbefore Read_q:\n");
        PSE::Read_q(q,y,ny,z,nz,alpha,"uvwP_251");
        PSE::printVecView(q);
        PSE::printScalar(y,ny);
        //PSE::Read_q(q,y,ny,z,nz,alpha);
        PetscPrintf(PETSC_COMM_WORLD,"\noutput:\n");
        PSE::printScalar(&alpha);
        // set A,b 
        PSE::set_A_and_B(y,ny,z,nz,A,B,Re,rho,alpha,m,omega);
        PSE::set_Euler_Advance(hx,A,B);
        PSE::Init_Vec(b,dim);
        PSE::set_b(B,q,b); //B*q->b
        // set BCs
        PSE::set_BCs(A,b,ny,nz);
        // view A
        //PSE::printMatView(A);
        //PSE::printMatASCII(A,"printMatASCII_dense.txt",PETSC_VIEWER_ASCII_DENSE);
        //PSE::printMatASCII(A);
        // solve Ax=b
        PSE::set_Vec(b); // assemble b
        PSE::Ax_b(A,qp1,b,dim);
        // view solution
        PSE::printVecView(q);
        PSE::printVecView(qp1);
        PSE::printVecASCII(q  ,"printVecq.txt");
        PSE::printVecASCII(qp1,"printVecqp1.txt");
        //PSE::printVecView(qexact);
        // free memory
        ierr = MatDestroy(&A);CHKERRQ(ierr);
        ierr = MatDestroy(&B);CHKERRQ(ierr);
        ierr = VecDestroy(&b);CHKERRQ(ierr);
        ierr = VecDestroy(&q);CHKERRQ(ierr);
        ierr = VecDestroy(&qp1);CHKERRQ(ierr);
    }
    if(0){ // compare against matmult and Ax_b using data
        // init
        Mat A,B;
        Vec b,q,Aq;
        PetscScalar Re=6000.,rho=1.,alpha,m=1.,omega=0.27;
        PetscScalar none=-1.0;
        PetscReal norm;
        PetscInt ny=201,nz=6;
        PetscScalar y[ny],z[nz];
        PetscScalar hx=0.001;
        PetscInt dim=ny*nz*4;
        PSE::Init_Vec(q,dim);
        // read in q,y,z,alpha from binary files
        PSE::Read_q(q,y,ny,z,nz,alpha);
        PetscPrintf(PETSC_COMM_WORLD,"\noutput:\n");
        PSE::printScalar(&alpha);
        // set A,b 
        PSE::set_A_and_B(y,ny,z,nz,A,B,Re,rho,alpha,m,omega);
        //PSE::set_Euler_Advance(hx,A,B);
        PSE::Init_Vec(b,dim);
        PSE::set_b(B,q,b); //B*q->b
        // set BCs
        PSE::set_BCs(A,b,ny,nz);
        PSE::set_Vec(b); // assemble b
        // view A
        //PSE::printMatView(A);
        //PSE::printMatASCII(A,"printMatASCII_dense.txt",PETSC_VIEWER_ASCII_DENSE);
        //PSE::printMatASCII(A);
        // Check if Aq = 0 or not
        PSE::Init_Vec(Aq,dim);
        MatMult(A,q,Aq);
        // view solution
        PSE::printVecView(q);
        PSE::printVecView(Aq);
        PSE::printVecASCII(Aq,"Aq.txt");
        //PSE::printVecView(qexact);
        // free memory
        ierr = MatDestroy(&A);CHKERRQ(ierr);
        ierr = MatDestroy(&B);CHKERRQ(ierr);
        ierr = VecDestroy(&b);CHKERRQ(ierr);
        ierr = VecDestroy(&q);CHKERRQ(ierr);
        ierr = VecDestroy(&Aq);CHKERRQ(ierr);
    }

    ierr = PetscFinalize();
    return ierr;


}
