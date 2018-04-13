#ifndef INIT_X_H
#define INIT_X_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief Initialize a vector x to be of size n
     * \return 0 if successful
     */
    int Init_x(
            Vec &x,         ///< x vector to initialize in MPI
            PetscInt &n      ///< size of vector
            );
}

#endif
