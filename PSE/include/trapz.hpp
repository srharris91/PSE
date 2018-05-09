#ifndef TRAPZ_H
#define TRAPZ_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief trapezoidal rule on Vec with ny values
     * \return 0 if successful
     */
    int trapz(
            Vec &q,             ///< Vector to integrate
            const PetscInt &n,  ///< size of vector
            PetscScalar &I,    ///< output value of integration 
            const PetscScalar &Deltay=1. ///< Domain size (b-a) in \f$ I = \int_a^b f(x) dx\f$
            );
    /**
     * \brief trapezoidal rule on Vec with ny,nz values
     * \return 0 if successful
     */
    int trapz(
            Vec &q,                         ///< Vector to integrate
            const PetscInt &ny,             ///< size of vector in ny (length(q) = 4*ny*nz)
            const PetscInt &nz,             ///< size of vector in nz (length(q) = 4*ny*nz)
            PetscScalar &I,                 ///< output value of integration 
            const PetscScalar &Deltay=1.,   ///< Domain size (b-a) in y-direction
            const PetscScalar &Deltaz=1.    ///< Domain size (b-z) in z-direction
            );



}
#endif
