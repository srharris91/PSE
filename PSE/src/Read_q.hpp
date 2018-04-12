#ifndef READ_Q_H
#define READ_Q_H
#include <petscksp.h>

namespace PSE
{
    /** 
     * Read q from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            PetscScalar *RHS_True,  ///< PetscScalar vector of true RHS values in Ax=RHS matrix
            int n                   ///< int size of vector
            );
}
#endif
