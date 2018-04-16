#ifndef PRINT_H
#define PRINT_H

#include <petscksp.h>
namespace PSE
{
    /** \brief Print PetscScalar variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printScalar(
            PetscScalar x[],    ///< PetscScalar array to print to screen
            PetscInt n=1,              ///< size of scalar array to print
            char const name[]="x",    ///< name of variable to output default to 'x'
            PetscViewer viewer=PETSC_VIEWER_STDOUT_WORLD ///< format for viewer
            ); 

    /** \brief Print Vec from PETSc type variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printVec(
            Vec &x,    ///< Vec array to print to screen
            PetscInt n,              ///< size of scalar array to print
            char const name[]="x"     ///< name of variable to output default to 'x'
            ); 
    /** \brief Print PetscInt variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printInt(
            PetscInt x[],    ///< PetscScalar array to print to screen
            PetscInt n=1,              ///< size of scalar array to print
            char const name[]="x"     ///< name of variable to output default to 'x'
            ); 
    /** \brief view PetscVec variable to screen
     *   Open an X-window viewer.  Note that we specify the same communicator
     *    for the viewer as we used for the distributed vector (PETSC_COMM_WORLD).
     *
     *    - Helpful runtime option:
     *
     *    -draw_pause <pause> : sets time (in seconds) that the
     *    program pauses after PetscDrawPause() has been called
     *    (0 is default, -1 implies until user input).
     */
    void printVecView( 
            Vec &x,             ///< PetscScalar array to print to screen
            PetscInt n,         ///< size of scalar array to print
            char const name[]="x viewer",  ///< name of variable to output default to 'x'
            PetscViewerFormat format=PETSC_VIEWER_DEFAULT ///< format for viewer
            );
    /** \brief view Petsc Mat variable to screen
     *   Open an X-window viewer.  Note that we specify the same communicator
     *    for the viewer as we used for the distributed vector (PETSC_COMM_WORLD).
     *
     *    - Helpful runtime option:
     *
     *    -draw_pause <pause> : sets time (in seconds) that the
     *    program pauses after PetscDrawPause() has been called
     *    (0 is default, -1 implies until user input).
     */
    void printMatView( 
            Mat &A,             ///< Petsc Mat nxn matrix to print to screen
            PetscInt n,         ///< size of nxn matrix to print
            char const name[]="A viewer",  ///< name of variable to output default to 'A'
            PetscViewerFormat format=PETSC_VIEWER_DEFAULT ///< format for viewer
            );
}

#endif
