#ifndef SET_BCS_H
#define SET_BCS_H
#include <petscksp.h>
namespace PSE
{
    /** 
     * 
     * \brief sets Boundary Condition equations in Matrix A and vector b in Ax=b linear system for solving the PSE equations
     * \return 0 if successful
     */
    PetscInt set_BCs(
            Mat &A,                 ///< [in,out] A matrix (already initialized)
            Vec &b,                 ///< [in,out] b Vec (already initialized)
            const PetscInt &ny,     ///< [in] size of y array
            const PetscInt &nz     ///< [in] size of z array
            );
    /** 
     * 
     * \brief sets Boundary Condition equations in Matrix A and Matrix B in A+ B*dqdx=0 linear system for solving the PSE equations
     * \return 0 if successful
     */
    PetscInt set_BCs(
            Mat &A,                 ///< [in,out] A matrix (already initialized)
            Mat &B,                 ///< [in,out] B Matrix (already initialized)
            const PetscInt &ny,     ///< [in] size of y array
            const PetscInt &nz     ///< [in] size of z array
            );
}

#endif

