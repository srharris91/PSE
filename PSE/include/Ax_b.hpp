#ifndef AX_B_H
#define AX_B_H
#include <petscksp.h>

namespace PSE
{
    /** 
     * \brief Solves \f$x\f$ in \f$Ax=b\f$ linalg problem using PETSc in parallel
     * \return 0 if successful
     */
    PetscInt Ax_b(
            const PetscScalar* const* Ain,  ///< input of 2D \f$A\f$ matrix in \f$Ax=b\f$ (must be dynamic pointer to pointer)
            Vec &x,             ///< vector of \f$x\f$ in \f$Ax=b\f$
            const PetscScalar bin[],  ///< vector of \f$b\f$ in \f$Ax=b\f$
            const PetscInt &n          ///< n size of vectors or nxn matrix \f$A\f$
            );
    /** 
     * \brief Solves \f$x\f$ in \f$Ax=b\f$ linalg problem using PETSc in parallel
     * \return 0 if successful
     */
    PetscInt Ax_b(
            const Mat &A,              ///< input of A matrix in \f$Ax=b\f$
            Vec &x,             ///< vector of \f$x\f$ in \f$Ax=b\f$
            const Vec &b,              ///< vector of \f$b\f$ in \f$Ax=b\f$
            const PetscInt &n          ///< n size of vectors or nxn matrix \f$A\f$
            );
}
#endif
