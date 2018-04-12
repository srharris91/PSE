#ifndef BASE_FLOW_H
#define BASE_FLOW_H
#include <petscksp.h>
/** 
 * \brief set plane channel flow velocity variables U, U' and U''
 * \usage use this once to set the vectors with the velocity of plane channel flow (max u=1)
 * \return 0 if successful
 */
int base_flow(
        PetscScalar U[],    //!< velocity along channel axis
        PetscScalar Uy[],   //!< first derivative of velocity
        PetscScalar Uyy,    //!< second derivative of velocity 
        PetscScalar y[],    //!< input y vector (usually linspace(-1,1,n)
        int n,              //!< size of vectors
        PetscBool output_full=PETSC_FALSE, //!< do we output the full vector, or just the inner vector (no the edges)
        PetscInt Dim=2      //!< Dimension of array to be put return (may not need this now
        );
#endif
