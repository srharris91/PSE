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
            const PetscScalar &hx,  ///< delta x distance
            Mat &A,                 ///< A matrix set from set_A_and_B
            Mat &B                  ///< B matrix set from set_A_and_B
            );
}

#endif
