#ifndef INIT_VEC_H
#define INIT_VEC_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief Initialize a vector x to be of size n
     * \return 0 if successful
     */
    int Init_Vec(
            Vec &x,         ///< x vector to initialize in MPI
            const PetscInt &n      ///< size of vector
            );
}

#endif
