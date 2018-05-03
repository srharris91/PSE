#ifndef TRAPZ_H
#define TRAPZ_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief trapezoidal rule on Vec
     * \return 0 if successful
     */
    int trapz(
            Vec &q,             ///< Vector to integrate
            const PetscInt &n,  ///< size of vector
            PetscScalar &I,    ///< output value of integration 
            const PetscScalar &Deltax=1. ///< Domain size (b-a) in \f[ I = \int_a^b f(x) dx\f]
            );



}
#endif
