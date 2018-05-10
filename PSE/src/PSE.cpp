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
    if(0){ // advance q one step, and check growth
        // init
        Mat A,B;
        Vec b,q,qp1;
        PetscScalar Re=6000.,rho=1.,alpha,m=1.,omega=0.27;
        PetscInt ny=101,nz=6;
        PetscScalar y[ny],z[nz];
        PetscScalar hx=2.5;
        PetscInt dim=ny*nz*4;
        PSE::Init_Vec(q,dim);
        // read in q,y,z,alpha from binary files
        PSE::Read_q(q,y,ny,z,nz,alpha,"../OrrSommerfeld_and_primitive/uvwP_101");
        //PSE::printVecView(q);
        //PSE::printScalar(y,ny);
        //PSE::Read_q(q,y,ny,z,nz,alpha);
        PetscPrintf(PETSC_COMM_WORLD,"\noutput:\n");
        PSE::printScalar(&alpha);
        // set A,b 
        PSE::set_A_and_B(y,ny,z,nz,A,B,Re,rho,alpha,m,omega);
        // set BCs in A and B
        PSE::set_BCs(A,B,ny,nz);
        // set up Euler advancing matrices
        PSE::set_Euler_Advance(hx,A,B);
        // set b vector from b=B*q
        PSE::Init_Vec(b,dim);
        PSE::set_b(B,q,b); //B*q->b
        // set BCs alternate BC setting
        //PSE::set_BCs(A,b,ny,nz);
        // view A
        //PSE::printMatView(A);
        //PSE::printMatASCII(A,"printMatASCII_dense.txt",PETSC_VIEWER_ASCII_DENSE);
        //PSE::printMatASCII(A);
        // solve Ax=b
        PSE::set_Vec(b); // assemble b
        PSE::Ax_b(A,qp1,b,dim);
        // view solution
        //PSE::printVecView(q);
        //PSE::printVecView(qp1);
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
    if(0){ // compare against matmult and Ax_b using data (should be exactly equal to OSS equations from python scripts)
        // init
        Mat A,B;
        Vec b,q,Aq;
        PetscScalar Re=6000.,rho=1.,alpha,m=1.,omega=0.27;
        PetscInt ny=101,nz=6;
        PetscScalar y[ny],z[nz];
        PetscInt dim=ny*nz*4;
        PSE::Init_Vec(q,dim);
        // read in q,y,z,alpha from binary files
        //PSE::Read_q(q,y,ny,z,nz,alpha);
        PSE::Read_q(q,y,ny,z,nz,alpha,"../OrrSommerfeld_and_primitive/uvwP_101");
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
        PSE::printMatASCII(A,"A.m",PETSC_VIEWER_ASCII_MATLAB);
        // Check if Aq = 0 or not
        PSE::Init_Vec(Aq,dim);
        MatMult(A,q,Aq);
        // view solution
        //PSE::printVecView(q);
        //PSE::printVecView(Aq);
        PSE::printVecASCII(Aq,"Aq.txt");
        // free memory
        ierr = MatDestroy(&A);CHKERRQ(ierr);
        ierr = MatDestroy(&B);CHKERRQ(ierr);
        ierr = VecDestroy(&b);CHKERRQ(ierr);
        ierr = VecDestroy(&q);CHKERRQ(ierr);
        ierr = VecDestroy(&Aq);CHKERRQ(ierr);
    }
    if(0){ // get subvector test and trapz
        //PetscScalar s[]={-3,-2,-1,0};
        PetscInt n=202;
        PetscScalar *s=new PetscScalar[n];
        for(int i=0; i<n; ++i) s[i]=i;
        Vec sVec,sVecsub;
        PSE::Init_Vec(sVec,n);
        PSE::set_Vec(s,n,sVec);
        PSE::set_Vec(sVec);
        PSE::printVecView(sVec);

        // set subvector
        PSE::set_Vec(sVec,0,101,sVecsub);

        // trapz
        PetscScalar trap_value;
        PSE::trapz(sVecsub,101,trap_value,100); // 5000 as expected! 

        // view
        PSE::printVecView(sVecsub);
        PSE::printScalar(&trap_value); // x=5000 as expected!

        // destroy
        //VecRestoreSubVector(sVec,is,&sVecsub);
        ierr = VecDestroy(&sVec);CHKERRQ(ierr);
        ierr = VecDestroy(&sVecsub);CHKERRQ(ierr);
        delete[] s;
        //delete[] idx;
        //ierr = ISDestroy(&is);CHKERRQ(ierr);

    }
    if(1){ // advance q many steps, and check growth, and iterate on alpha
        // init
        Mat A,B;
        Vec b,q,qp1;
        //PetscScalar Re=6000.,rho=1.,alpha,m=1.,omega=0.27;
        PetscScalar Re=2000.,rho=1.,alpha,m=1.,omega=0.3;
        PetscInt ny=101,nz=6;
        PetscScalar y[ny],z[nz];
        PetscScalar hx=2.5;
        PetscInt dim=ny*nz*4;
        PSE::Init_Vec(q,dim);
        PSE::Init_Vec(qp1,dim);
        // read in and set matrices
        PSE::Read_q(q,y,ny,z,nz,alpha,"../OrrSommerfeld_and_primitive/uvwP_101_stable");// read in q,y,z,alpha from binary files
        VecCopy(q,qp1); // qp1=q for initial guess

        //VecSetValue(q,ny+6,0.00001,ADD_VALUES);// add a little v disturbance
        PSE::set_Vec(q); // assemble again
        //PSE::printVecView(q);
        char filename[100];
        sprintf(filename,"printVecq0.txt");
        PSE::printVecASCII(q  ,filename);
        PSE::printScalar(&alpha,1,"original alpha");
        PetscScalar Ialpha = alpha*hx;
        for(int i=1; i<15; ++i){
            PSE::set_A_and_B(y,ny,z,nz,A,B,Re,rho,alpha,m,omega);// set A,b 
            PSE::set_BCs(A,B,ny,nz);        // set BCs in A and B
            PSE::set_Euler_Advance(hx,A,B); // set up Euler advancing matrices
            PSE::Init_Vec(b,dim);           // set b vector from b=B*q
            PSE::set_b(B,q,b);              //B*q->b
            PSE::set_Vec(b); // assemble b
            // solve Ax=b
            PSE::Ax_b(A,qp1,b,dim);
            // view solution
            //PSE::printVecView(q);
            //PSE::printVecView(qp1);
            // closure?
            PSE::update_Closure(q,qp1,ny,nz,hx,Ialpha,alpha);
            //PSE::printVecView(qp1);
            VecCopy(qp1,q);
            PSE::set_Vec(q); //assemble
            sprintf(filename,"printVecq%d.txt",i);
            PSE::printVecASCII(q  ,filename);
            PetscPrintf(PETSC_COMM_WORLD,"Marched %d\n",i);
        }
        //PSE::printVecASCII(q  ,"printVecq.txt");
        //PSE::printVecASCII(qp1,"printVecqp1.txt");
        // free memory
        ierr = MatDestroy(&A);CHKERRQ(ierr);
        ierr = MatDestroy(&B);CHKERRQ(ierr);
        ierr = VecDestroy(&b);CHKERRQ(ierr);
        ierr = VecDestroy(&q);CHKERRQ(ierr);
        ierr = VecDestroy(&qp1);CHKERRQ(ierr);
    }
    if(0){ // read in variables and update_Closure step a lot
        // note! Must force many iterations in update_Closure.cpp file for this to print and plot
        // init
        Vec q,qp1;
        //PetscScalar Re=6000.,rho=1.,alpha,m=1.,omega=0.27;
        PetscScalar Re=2000.,rho=1.,alpha,m=1.,omega=0.3;
        PetscInt ny=101,nz=6;
        PetscScalar y[ny],z[nz];
        PetscScalar hx=2.5;
        PetscInt dim=ny*nz*4;
        PSE::Init_Vec(q,dim);
        PSE::Init_Vec(qp1,dim);
        // read in and set matrices
        PSE::Read_q(q,y,ny,z,nz,alpha,"../OrrSommerfeld_and_primitive/uvwP_101_stable");// read in q,y,z,alpha from binary files
        //VecSetValue(q,ny+6,0.00001,ADD_VALUES);// add a little v disturbance
        PSE::set_Vec(q); // assemble again
        //PSE::printVecView(q);
        char filename[100];
        sprintf(filename,"printVecq0.txt");
        PSE::printVecASCII(q  ,filename);
        PSE::printScalar(&alpha,1,"original alpha");
        PetscScalar Ialpha = alpha*hx;

        // closure updates
        VecCopy(q,qp1); // qp1=q
        PSE::set_Vec(qp1);
        //VecSetValue(qp1,6,1,ADD_VALUES);// add a little v disturbance
        VecScale(qp1,1.+1e-1);
        //VecSetValue(qp1,7,-1e-10,ADD_VALUES);// add a little v disturbance
        PSE::set_Vec(qp1);
        PSE::update_Closure(q,qp1,ny,nz,hx,Ialpha,alpha);
        //PSE::printVecView(qp1);
        sprintf(filename,"printVecq%d.txt",1);
        PSE::printVecASCII(qp1  ,filename);

        // free memory
        ierr = VecDestroy(&q);CHKERRQ(ierr);
        ierr = VecDestroy(&qp1);CHKERRQ(ierr);
    }


    ierr = PetscFinalize();
    return ierr;


}
