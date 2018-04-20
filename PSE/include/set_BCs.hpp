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
            Mat &A,                 ///< A matrix
            Vec &b,                 ///< b Vec
            const PetscInt &ny,     ///< size of y array
            const PetscInt &nz     ///< size of z array
            );
}

#endif

