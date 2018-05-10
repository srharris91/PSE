#ifndef SET_EULER_ADVANCE_H
#define SET_EULER_ADVANCE_H
#include <petscksp.h>
namespace PSE
{
    /**
     * \brief matrix multiplication for Euler stepping, A+= 1/hx * B and B*=1/hx
     * \return 0 if successful
     *
     */
    PetscInt set_Euler_Advance(
            const PetscScalar &hx,  ///< [in] delta x distance
            Mat &A,                 ///< [in,out] A matrix set from set_A_and_B (already initialized)
            Mat &B                  ///< [in,out] B matrix set from set_A_and_B (already initialized)
            );
}

#endif
