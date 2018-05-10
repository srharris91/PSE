#ifndef BASE_FLOW_H
#define BASE_FLOW_H
#include <petscksp.h>

namespace PSE
{
    /** 
     * \brief set plane channel flow velocity variables U, U' and U''
     * \usage use this once to set the vectors with the velocity of plane channel flow (max u=1)
     * \return 0 if successful
     */
    int base_flow(
            PetscScalar U[],        ///< [out] velocity along channel axis
            PetscScalar Uy[],       ///< [out] first derivative of velocity
            const PetscScalar y[],  ///< [in] input y vector (usually linspace(-1,1,n)
            const PetscInt &n       ///< [in] size of vectors
            );
    /** 
     * \brief set plane channel flow velocity variables U, U' and U''
     * \usage use this once to set the vectors with the velocity of plane channel flow (max u=1)
     * \return 0 if successful
     */
    int base_flow(
            Mat &U,                 ///< [out] velocity along channel axis (uninitialized)
            Mat &Uy,                ///< [out] first derivative of velocity (uninitialized)
            const PetscScalar y[],  ///< [in] input y vector (usually linspace(-1,1,n)
            const PetscInt &n       ///< [in] size of vectors
            );
            }
#endif
