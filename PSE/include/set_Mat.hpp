#ifndef SET_MAT_H
#define SET_MAT_H
#include <petscksp.h>
namespace PSE{
    /**
     * \brief set a matrix from PetscScalar 2D matrix to PETSc Mat type in parallel
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar* const* Ain,  //< Matrix set on all processors as 2D array
            const PetscInt &n,        //< size of nxn matrix
            Mat &A              //< Mat to SetValues and output in parallel
            );

    /**
     * \brief set a diagonal of a matrix from scalar
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &diag,    //< scalar value to set as diagonal in matrix
            const PetscInt &n,          //< size of nxn matrix
            Mat &A,                     //< Mat to SetValues and output in parallel
            const PetscInt &k=0         //< diagonal offset in matrix (k=0 is main diagonal, k>0 is above main diagonal, k<0 is below main diagonal) default is 0
            );
}


#endif
