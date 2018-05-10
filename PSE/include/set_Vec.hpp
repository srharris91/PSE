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
            const PetscScalar *bin,   ///< [in] array set on all processors as 1D array
            const PetscInt &n,        ///< [in] size of n array
            Vec &b,             ///< [in,out] Vec to SetValues and output in parallel (already initialized)
            const PetscBool &parallel=PETSC_TRUE ///< [in] set in parallel
            );

    /**
     * \brief Assemble b vector
     * \return 0 if successful
     *
     */
    PetscInt set_Vec(
            Vec &b ///< [in,out] array to assemble on all processors
            );
    /**
     * \brief set a vector from PetscScalar to PETSc Vec in location n
     * \return 0 if successful
     *
     */
    PetscInt set_Vec(
            const PetscScalar &bin,   ///< [in] scalar to set into Vec
            const PetscInt &n,        ///< [in] location to put into Vec
            Vec &b                    ///< [in,out] Vec to SetValues (already initialized)
            );
    /**
     * \brief set a subvector from larger vector from low to hi indices
     * \return 0 if successful
     *
     */
    PetscInt set_Vec(
            const Vec &inVec,           ///< [in] larger vec to copy values from
            const PetscInt &low,        ///< [in] lower bound
            const PetscInt &hi,        ///< [in] upper bound
            Vec &b                    ///< [in,out] Vec to SetValues and output (already initialized)
            );
}


#endif
