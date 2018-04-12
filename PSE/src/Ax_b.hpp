#ifndef AX_B_H
#define AX_B_H
#include <petscksp.h>

namespace PSE
{
    /** 
     * \brief Solves \f$x\f$ in \f$Ax=b\f$ linalg problem using PETSc in parallel
     * \return ierr from Petsc calls 
     */
    PetscErrorCode Ax_b(
            PetscScalar Ain1[], ///< first row input of 2D \f$A\f$ matrix in \f$Ax=b\f$
            PetscScalar xout[], ///< vector of \f$x\f$ in \f$Ax=b\f$
            PetscScalar bin[],  ///< vector of \f$b\f$ in \f$Ax=b\f$
            PetscInt n,         ///< n size of vectors or nxn matrix \f$A\f$
            int argc,           ///< argc from int main()
            char **args         ///< args from int main
            );
}
#endif
