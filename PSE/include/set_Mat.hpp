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
            const PetscScalar* const* Ain,  ///< [in] Matrix set on all processors as 2D array
            const PetscInt &n,              ///< [in] size of nxn matrix
            Mat &A,                         ///< [out] Mat to SetValues and output in parallel (uninitialized)
            const InsertMode &addv=ADD_VALUES ///< [in] insert values or add values to matrix A
            );

    /**
     * \brief set a diagonal of a matrix from scalar
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &diag,    ///< [in] scalar value to set as diagonal in matrix
            const PetscInt &n,          ///< [in] size of nxn matrix
            Mat &A,                     ///< [in,out] Mat to SetValues and output in parallel (already initialized)
            const PetscInt &k=0,        ///< [in] diagonal offset in matrix (k=0 is main diagonal, k>0 is above main diagonal, k<0 is below main diagonal) default is 0
            const PetscBool &parallel=PETSC_TRUE, ///< [in] set the matrix using Istart and Iend
            const InsertMode &addv=ADD_VALUES ///< [in] insert values or add values to matrix A
            );
    /**
     * \brief Assemble matrix
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            Mat &A  ///< [in,out] Mat to assemble (already initialized)
            );

    /**
     * \brief set a single value in matrix at row,col
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &value,   ///< [in] value to set in Mat A
            const PetscInt &row,        ///< [in] specified global row in Mat A
            const PetscInt &col,        ///< [in] specified global col in Mat A
            Mat &A,                     ///< [in,out] Matrix A to set value (already initialized)
            const PetscBool &parallel=PETSC_TRUE,  ///< [in] do this in parallel?
            const InsertMode &addv=ADD_VALUES ///< [in] insert values or add values to matrix A
            );

    /**
     * \brief set a submatrix Asub in a matrix A
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &a,       ///< [in] premultple to Asub
            const Mat &Asub,            ///< [in] sub matrix to set into mat
            const PetscInt &nsub,       ///< [in] nxn size of square sub matrix Asub
            Mat &A,                     ///< [in,out] Matrix A to set value (already initialized)
            const PetscInt &n,          ///< [in] size of square Matrix A
            const PetscInt &rowoffset=0,  ///< [in] start inserting Asub in A starting at rowoffset row
            const PetscInt &coloffset=0,  ///< [in] start inserting Asub in A starting at coloffset column
            const InsertMode &addv=ADD_VALUES ///< [in] insert values or add values to matrix A
            );

    /**
     * \brief set a row of PetscScalar to a matrix to PETSc Mat 
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar* Ain,         ///< [in] array to set in row of matrix
            const PetscInt &row,            ///< [in] row to set in matrix
            const PetscInt &ncols,          ///< [in] ncols of array
            const PetscInt cols[],          ///< [in] cols to set in matrix
            Mat &A,                         ///< [in,out] Mat to SetValues and output in parallel (already initialized)
            const InsertMode &addv=ADD_VALUES ///< [in] insert values or add values to matrix A
            );

    /**
     * \brief set a Dz submatrix into A
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &a,       ///< [in] premultple to Dz
            const Mat &Dz,            ///< [in] sub matrix to set into mat in z direction
            const PetscInt &nz,       ///< [in] nzxnz size of square sub matrix Dz
            const PetscInt &ny,          ///< [in] length of y array
            const PetscInt &zi,         ///< [in] which z plane you are looking at
            Mat &A,                     ///< [in,out] Matrix A to set value (already initialized)
            const PetscInt &n,          ///< [in] size of square Matrix A
            const PetscInt &rowoffset=0,  ///< [in] start inserting Asub in A starting at rowoffset row
            const PetscInt &coloffset=0,  ///< [in] start inserting Asub in A starting at coloffset column
            const InsertMode &addv=ADD_VALUES ///< [in] insert values or add values to matrix A
            );
}


#endif
