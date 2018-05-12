#ifndef SET_A_AND_B_H
#define SET_A_AND_B_H
#include <petscksp.h>
#include "class.hpp"
namespace PSE
{
    /** 
     * 
     * \brief set A and B matrix for PSE equations
     * \return 0 if successful
     */
    PetscInt set_A_and_B(
            PSE &data   ///< data class containing all matrices, vectors, scalars, and problem size
            );
}

#endif

