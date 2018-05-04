#ifndef UPDATE_CLOSURE_H
#define UPDATE_CLOSURE_H
#include <petscksp.h>
namespace PSE{
    /** 
     * \brief calculate the closure equation \f[ \int_\Omega \hat\textbf{q}^\dagger \hat\textbf{q}_x dy \f] and output the result
     * \return 0 if successful
     */
    int update_Closure(
            const Vec &q,                       ///< input array from previous iteration \f[ \hat\textbf{q}_i \f]
            Vec &qp1,                           ///< input array from previous iteration \f[ \hat\textbf{q}_{i+1} \f]
            const PetscInt &ny,                 ///< size of arrays in y dimension
            const PetscInt &nz,                 ///< size of arrays in z dimension
            const PetscScalar &dx,                 ///< x marching step distance
            PetscScalar &Ialpha,                ///< integral value from previous step \f[ \int_0^{x_i} \alpha(\xi) d\xi \f]
            PetscScalar &alpha,                 ///< initial alpha value of eigenfunction (constant from last step initially)
            const PetscScalar &tol=1e-16+1e-18*PETSC_i, ///< tolerance on norm equation
            const PetscScalar &Deltax=2         ///< height of channel (default to 2)
            );

}

#endif
