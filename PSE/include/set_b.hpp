#ifndef SET_b_H
#define SET_b_H
#include <petscksp.h>
namespace PSE
{
    /** 
     * 
     * \brief set b vector from B and q as b=B*q
     * \return 0 if successful
     */
    PetscInt set_b(
            const Mat &B,                 ///< B matrix 
            const Vec &qn,          ///< qn vector to multiply
            Vec &b                 ///< b Vec (initialized)
            );
}

#endif

