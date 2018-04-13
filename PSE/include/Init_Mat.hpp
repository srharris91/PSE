#ifndef INIT_MAT_H
#define INIT_MAT_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief Initialize a Matrix A to be of size nxn
     * \return 0 if successful
     */
    int Init_Mat(
            Mat &A,         ///< A matrix to initialize in MPI
            const PetscInt &n      ///< size of vector
            );
}

#endif
