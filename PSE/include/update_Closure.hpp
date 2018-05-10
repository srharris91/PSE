#ifndef UPDATE_CLOSURE_H
#define UPDATE_CLOSURE_H
#include <petscksp.h>
namespace PSE{
    /** 
     * \brief calculate the closure equation \f$ \int_\Omega \hat{\textbf{q}}^\dagger \hat{\textbf{q}}_x dy \f$ and output the result
     * This is part of the solution procedure shown in the bottom half of this diagram
     * \image html Diagram_closure.png "Solution Procedure" width=10cm
     *
     *
     *
     * \return 0 if successful
     */
    int update_Closure(
            const Vec &q,                       ///< [in] array from previous iteration \f$ \hat{\textbf{q}}_i \f$
            Vec &qp1,                           ///< [in,out] array from current iteration \f$ \hat{\textbf{q}}_{i+1} \f$
            const PetscInt &ny,                 ///< [in] size of arrays in y dimension
            const PetscInt &nz,                 ///< [in] size of arrays in z dimension
            const PetscScalar &dx,              ///< [in] x marching step distance
            PetscScalar &Ialpha,                ///< [in,out] integral value from previous step \f$ \int_0^{x_i} \alpha(\xi) d\xi \f$
            PetscScalar &alpha,                 ///< [in,out] initial alpha value of eigenfunction (constant from last step initially)
            const PetscReal &tol=1e-7,         ///< [in] tolerance on norm equation
            const PetscScalar &Deltay=2,        ///< [in] height of channel in y (default to 2)
            const PetscScalar &Deltaz=1         ///< [in] width of channel in z(default to 1)
            );

}

#endif
