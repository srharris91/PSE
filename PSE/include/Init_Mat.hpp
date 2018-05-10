#ifndef INIT_MAT_H
#define INIT_MAT_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief Initialize a Matrix A to be of size nxn
     * \return 0 if successful
     */
    int Init_Mat(
            Mat &A,         ///< [out] A matrix to initialize in MPI (uninitialized)
            const PetscInt &n      ///< [in] number of global rows and columns
            );
    /**
     * \brief Initialize a Matrix A to be of size mxn
     * \return 0 if successful
     */
    int Init_Mat(
            Mat &A,         ///< [out] A matrix to initialize in MPI (uninitialized)
            const PetscInt &m,      ///< [in] number of global rows
            const PetscInt &n      ///< [in] number of global columns
            );
}

#endif
