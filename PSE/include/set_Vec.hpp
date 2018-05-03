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
            Vec &b,             ///< Vec to SetValues and output in parallel
            const PetscBool &parallel=PETSC_TRUE ///< set in parallel
            );

    /**
     * \brief Assemble b vector
     * \return 0 if successful
     *
     */
    PetscInt set_Vec(
            Vec &b ///< array to assemble on all processors
            );
    /**
     * \brief set a vector from PetscScalar to PETSc Vec in location n
     * \return 0 if successful
     *
     */
    PetscInt set_Vec(
            const PetscScalar &bin,   ///< scalar to set into Vec
            const PetscInt &n,        ///< location to put into Vec
            Vec &b                    ///< Vec to SetValues
            );
    /**
     * \brief set a subvector from larger vector from low to hi indices
     * \return 0 if successful
     *
     */
    PetscInt set_Vec(
            const Vec &inVec,   ///< larger vec to copy values from
            const PetscInt &low,        ///< lower bound
            const PetscInt &hi,        ///< upper bound
            Vec &b                    ///< Vec to SetValues and output
            );
}


#endif
