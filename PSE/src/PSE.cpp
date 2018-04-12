
static char help[] = "Solves a linear system in parallel with KSP.\n\n";

/*T
   Concepts: KSP^solving the PSE equations
   Concepts: complex numbers;
   Concepts: Parabolic Stability Equations for Transition modeling
   Processors: n
T*/

/*
   Description: Solves a complex linear system in parallel with KSP.

   The model problem:
      Solve PSE equation in the channel flow
      Dirichlet b.c.'s on all sides
      Use the 3-D, finite difference stencil.

   Compiling the code:
      This code uses the complex numbers version of PETSc, so configure
      must be run to enable this
*/

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include <petscksp.h>
#include <petscsys.h>
#include <stdio.h>
#include <iostream>

#include "PSE.hpp"

int main(int argc,char **args)
{
    // init some variables
  PetscInt       n = 6;
  PetscErrorCode ierr;
  PetscScalar    x[n],x2[n]; 
  PetscScalar    Atest[n][6] = {
      {1,1,-2,1,3,-4},
      {2,-1,1,2,1,-3},
      {1,3,-3,-1,2,1},
      {5,2,-1,-1,2,1},
      {-3,-1,2,3,1,3},
      {4,3,1,-6,-3,-2} };
  PetscScalar    btest[] = {7,20,-15,-3,16,-27};

  // initialize Petsc
  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

  // solve Ax=b problem
  ierr = Ax_b(Atest[0],x,btest,n,argc,args); CHKERRQ(ierr);
  // output solutions
  for(int i=0; i<n; i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  x[%D] = %g + %g i\n",i,(double)PetscRealPart(x[i]),(double)PetscImaginaryPart(x[i]));CHKERRQ(ierr);}

  // solve Ax=b problem
  ierr = Ax_b(Atest[0],x2,btest,n,argc,args); CHKERRQ(ierr);
  // output solutions
  for(int i=0; i<n; i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  x2[%D] = %g + %g i\n",i,(double)PetscRealPart(x2[i]),(double)PetscImaginaryPart(x[i]));CHKERRQ(ierr);}






  // finalize
  ierr = PetscFinalize();
  return ierr;
}
