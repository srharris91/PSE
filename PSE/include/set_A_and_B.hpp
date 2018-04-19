#ifndef SET_A_AND_B_H
#define SET_A_AND_B_H
#include <petscksp.h>
namespace PSE
{
    /** 
     * 
     * \brief set A and B matrix for PSE equations
     * \return 0 if successful
     */
    PetscInt set_A_and_B(
            const PetscScalar y[],  ///< y array
            const PetscInt &ny,     ///< size of y array
            const PetscScalar z[],  ///< z array
            const PetscInt &nz,     ///< size of z array
            Mat &A,                 ///< A matrix (uninitialized)
            Mat &B,                 ///< B matrix (uninitialized)
            const PetscInt &order=4,  ///< order of accuracy of derivatives in finite difference
            const PetscBool &reduce_wall_order=PETSC_TRUE ///< reduce order of derivatives at the wall if true
            );
}

#endif

