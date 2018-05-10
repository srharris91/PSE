#ifndef SET_A_AND_B_ZI_H
#define SET_A_AND_B_ZI_H
#include <petscksp.h>
namespace PSE
{
    /** 
     * 
     * \brief set A and B matrix for PSE equations for zi plane
     * \return 0 if successful
     */
    PetscInt set_A_and_B_zi(
            const PetscScalar y[],  ///< [in] y array
            const PetscInt &ny,     ///< [in] size of y array
            const PetscScalar z[],  ///< [in] z array
            const PetscInt &nz,     ///< [in] size of z array
            Mat &A,                 ///< [out] A matrix (uninitialized)
            Mat &B,                 ///< [out] B matrix (uninitialized)
            const PetscScalar &Re,  ///< [in] Reynolds number
            const PetscScalar &rho,  ///< [in] rho density
            const PetscScalar &alpha,///< [in] alpha
            const PetscScalar &m,   ///< [in] m for omega
            const PetscScalar &omega,///< [in] omega value
            const Mat &Dy,           ///< [in] Dy matrix
            const Mat &Dyy,           ///< [in] Dyy matrix
            const Mat &Dz,           ///< [in] Dz matrix
            const Mat &Dzz,           ///< [in] Dzz matrix
            const Mat &I,           ///< [in] I identity matrix
            const Mat &U,           ///< [in] U base flow matrix
            const Mat &Uy,           ///< [in] Uy base flow matrix
            const PetscInt &zi=0,     ///< [in] zi plane of matrix
            const PetscInt &order=4,  ///< [in] order of accuracy of derivatives in finite difference
            const PetscBool &reduce_wall_order=PETSC_TRUE ///< [in] reduce order of derivatives at the wall if true
            );
}

#endif

