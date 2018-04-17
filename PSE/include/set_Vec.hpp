#ifndef SET_VEC_H
#define SET_VEC_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief set a vector from PetscScalar 1D vector to PETSc Vec type in parallel
     * \return 0 if successful
     *
     */
    PetscInt set_Vec(
            const PetscScalar *bin,   ///< array set on all processors as 1D array
            const PetscInt &n,        ///< size of n array
            Vec &b              ///< Vec to SetValues and output in parallel
            );
}


#endif
