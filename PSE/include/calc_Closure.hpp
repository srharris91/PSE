#ifndef CALC_CLOSURE_H
#define CALC_CLOSURE_H
#include <petscksp.h>
namespace PSE{
    /** 
     * \brief calculate the closure equation \f$ \int_{\Omega} \hat{\textbf{q}}^\dagger \hat{\textbf{q}}_x dy \f$ and output the result
     * \return 0 if successful
     */
    int calc_Closure(
            const Vec &q,                       ///< input array from previous iteration \f$ \hat{\textbf{q}}_i \f$
            Vec &qp1,                           ///< input array from previous iteration \f$ \hat{\textbf{q}}_{i+1} \f$
            const PetscInt &ny,                 ///< size of arrays in y dimension
            const PetscInt &nz,                 ///< size of arrays in z dimension
            PetscScalar &I,                     ///< output integral value
            const PetscScalar &Deltax=2         ///< height of channel (default to 2)
            );

}

#endif
