#ifndef CALC_CLOSURE_H
#define CALC_CLOSURE_H
#include <petscksp.h>
namespace PSE{
    /** 
     * \brief calculate the closure equation \f$ \int_{\Omega} \hat{\textbf{q}}^\dagger \hat{\textbf{q}}_x dy \f$ and output the result
     * \return 0 if successful
     */
    int calc_Closure(
            const Vec &q,                       ///< [in] input array from previous iteration \f$ \hat{\textbf{q}}_i \f$
            Vec &qp1,                           ///< [in,out] input array from previous iteration \f$ \hat{\textbf{q}}_{i+1} \f$ (already initialized)
            const PetscInt &ny,                 ///< [in] size of arrays in y dimension
            const PetscInt &nz,                 ///< [in] size of arrays in z dimension
            PetscScalar &I,                     ///< [out] output integral value
            const PetscScalar &Deltay=2,        ///< [in] height of channel in y(default to 2)
            const PetscScalar &Deltaz=1         ///< [in] width of channel in z(default to 1)
            );

}

#endif
