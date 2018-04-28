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
            const PetscScalar* const* Ain,  ///< Matrix set on all processors as 2D array
            const PetscInt &n,              ///< size of nxn matrix
            Mat &A,                         ///< Mat to SetValues and output in parallel
            const InsertMode &addv=ADD_VALUES ///< insert values or add values to matrix A
            );

    /**
     * \brief set a diagonal of a matrix from scalar
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &diag,    ///< scalar value to set as diagonal in matrix
            const PetscInt &n,          ///< size of nxn matrix
            Mat &A,                     ///< Mat to SetValues and output in parallel
            const PetscInt &k=0,        ///< diagonal offset in matrix (k=0 is main diagonal, k>0 is above main diagonal, k<0 is below main diagonal) default is 0
            const PetscBool &parallel=PETSC_TRUE, ///< set the matrix using Istart and Iend
            const InsertMode &addv=ADD_VALUES ///< insert values or add values to matrix A
            );
    /**
     * \brief Assemble matrix
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            Mat &A  ///< Mat to assemble
            );

    /**
     * \brief set a single value in matrix at row,col
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &value,   ///< value to set in Mat A
            const PetscInt &row,        ///< specified global row in Mat A
            const PetscInt &col,        ///< specified global col in Mat A
            Mat &A,                     ///< Matrix A to set value
            const PetscBool &parallel=PETSC_TRUE,  ///< do this in parallel?
            const InsertMode &addv=ADD_VALUES ///< insert values or add values to matrix A
            );

    /**
     * \brief set a submatrix Asub in a matrix A
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &a,       ///< premultple to Asub
            const Mat &Asub,            ///< sub matrix to set into mat
            const PetscInt &nsub,       ///< nxn size of square sub matrix Asub
            Mat &A,                     ///< Matrix A to set value
            const PetscInt &n,          ///< size of square Matrix A
            const PetscInt &rowoffset=0,  ///< start inserting Asub in A starting at rowoffset row
            const PetscInt &coloffset=0,  ///< start inserting Asub in A starting at coloffset column
            const InsertMode &addv=ADD_VALUES ///< insert values or add values to matrix A
            );

    /**
     * \brief set a row of PetscScalar to a matrix to PETSc Mat 
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar* Ain,         ///< array to set in row of matrix
            const PetscInt &row,            ///< row to set in matrix
            const PetscInt &ncols,          ///< ncols of array
            const PetscInt cols[],          ///< cols to set in matrix
            Mat &A,                         ///< Mat to SetValues and output in parallel
            const InsertMode &addv=ADD_VALUES ///< insert values or add values to matrix A
            );

    /**
     * \brief set a Dz submatrix into A
     * \return 0 if successful
     *
     */
    PetscInt set_Mat(
            const PetscScalar &a,       ///< premultple to Dz
            const Mat &Dz,            ///< sub matrix to set into mat in z direction
            const PetscInt &nz,       ///< nzxnz size of square sub matrix Dz
            const PetscInt &ny,          ///< length of y array
            const PetscInt &zi,         ///< which z plane you are looking at
            Mat &A,                     ///< Matrix A to set value
            const PetscInt &n,          ///< size of square Matrix A
            const PetscInt &rowoffset=0,  ///< start inserting Asub in A starting at rowoffset row
            const PetscInt &coloffset=0,  ///< start inserting Asub in A starting at coloffset column
            const InsertMode &addv=ADD_VALUES ///< insert values or add values to matrix A
            );
}


#endif
