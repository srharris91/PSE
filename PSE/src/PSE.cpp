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
    // test Ax=b solver twice
    if (1){
        // init some variables
        PetscInt       n = 6;
        PetscErrorCode ierr;
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
        PetscScalar    btest[] = {7,20,-15,-3,16,-27};

        // initialize Petsc
        ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
        PSE::Init_x(x,n);
        PSE::Init_x(x2,n);

        // solve Ax=b problem
        ierr = PSE::Ax_b(Atest2,x,btest,n); CHKERRQ(ierr);
        // output solutions
        PSE::printVec(x,n);

        // solve Ax=b problem
        ierr = PSE::Ax_b(Atest2,x2,btest,n); CHKERRQ(ierr);
        // output solutions
        PSE::printVec(x2,n,"x2");

        // finalize
        for (int i=0; i<n; i++) delete[] Atest2[i];
        delete[] Atest2;
        //ierr = PetscFinalize();
        //return ierr;
    }
    // test Read_q
    if (0){
        PetscInt n=6;
        PetscScalar q[n];
        PetscErrorCode ierr;
        ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
        PSE::Read_q(q,n);
        PSE::printScalar(q,n,"After reading vector = ");
        ierr = PetscFinalize();
        return ierr;
    }
    // test get_D_Coeffs
    if (1){
        Vec x;
        PetscScalar s[]={-3,-2,-1,0};
        PetscInt n=4;
        //PetscScalar output[n];
        PetscErrorCode ierr;

        ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

        PSE::Init_x(x,n);


        PSE::get_D_Coeffs(s,n,x);
        PSE::printVec(x,n,"after");

        ierr = VecDestroy(&x);CHKERRQ(ierr);
        ierr = PetscFinalize();
        return 0;
    }


}
