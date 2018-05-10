#ifndef SET_D_H
#define SET_D_H
#include <petscksp.h>
namespace PSE
{
    /** 
     * 
     * \brief set D matrix operator for specified order and derivative
     * \return 0 if successful
     */
    PetscInt set_D(
            const PetscScalar y[],  ///< [in] array of y values of channel
            const PetscInt &n,      ///< [in] length of y values
            Mat &output,            ///< [out] matrix(n-2 by n) dth derivative of order O(h^order) assuming uniform y spacing (uninitialized)
            const PetscInt &order=2,  ///< [in] order of accuracy desired (assuming even e.g. 2,4,6,...)
            const PetscInt &d=2 ,     ///< [in] dth derivative
            const PetscBool &periodic=PETSC_FALSE,  ///< [in] periodic boundary
            const PetscBool &reduce_wall_order=PETSC_TRUE  ///< [in] reduce the order of accuracy at the wall?
            );
}

#endif

