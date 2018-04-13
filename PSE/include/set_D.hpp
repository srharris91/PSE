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
            const PetscScalar y[],  ///< array of y values of channel
            const PetscInt &n,      ///< length of y values
            Mat &output,            ///< output matrix(n-2 by n) dth derivative of order O(h^order) assuming uniform y spacing
            const PetscInt &order=2,  ///< order of accuracy desired (assuming even e.g. 2,4,6,...)
            const PetscInt &d=2 ,     ///< dth derivative
            const PetscBool &reduce_wall_order=PETSC_TRUE, ///< reduce the order of accuracy at the wall?
            const PetscBool &output_full=PETSC_FALSE, ///< output full matrix
            const PetscBool &periodic=PETSC_FALSE   ///< periodic boundary
            );
}

#endif

