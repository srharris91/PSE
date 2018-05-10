#ifndef PRINT_H
#define PRINT_H

#include <petscksp.h>
namespace PSE
{
    /** \brief Print PetscScalar variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printScalar(
            const PetscScalar x[],  ///< [in] PetscScalar array to print to screen
            const PetscInt n=1,     ///< [in] size of scalar array to print
            char const name[]="x",  ///< [in] name of variable to output default to 'x'
            const PetscViewer viewer=PETSC_VIEWER_STDOUT_WORLD ///< [in] format for viewer
            ); 

    /** \brief Print Vec from PETSc type variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printVec(
            const Vec &x,                 ///< [in] Vec array to print to screen
            const PetscInt n,             ///< [in] size of scalar array to print
            char const name[]="x"   ///< [in] name of variable to output default to 'x'
            ); 
    /** \brief Print PetscInt variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printInt(
            const PetscInt x[],           ///< [in] PetscScalar array to print to screen
            const PetscInt n=1,           ///< [in] size of scalar array to print
            char const name[]="x"   ///< [in] name of variable to output default to 'x'
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
            const Vec &x,                 ///< [in] PetscScalar array to print to screen
            char const name[]="x viewer",  ///< [in] name of variable to output default to 'x'
            const PetscViewerFormat format=PETSC_VIEWER_DEFAULT ///< [in] format for viewer
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
            const Mat &A,           ///< [in] Petsc Mat nxn matrix to print to screen
            char const name[]="A viewer",  ///< [in] name of variable to output default to 'A'
            const PetscViewerFormat format=PETSC_VIEWER_DEFAULT ///< [in] format for viewer
            );
    /** \brief view Petsc Mat variable to ASCII file
     *    - Helpful runtime option:
     *
     *    -draw_pause <pause> : sets time (in seconds) that the
     *    program pauses after PetscDrawPause() has been called
     *    (0 is default, -1 implies until user input).
     */
    void printMatASCII( 
            const Mat &A,           ///< [in] Petsc Mat nxn matrix to print to screen
            char const name[]="printMatASCII.txt", ///< [in] filename to write to
            const PetscViewerFormat format=PETSC_VIEWER_DEFAULT ///< [in] format for viewer
            );
    /** \brief view Petsc Vec variable to ASCII file
     *    - Helpful runtime option:
     *
     *    -draw_pause <pause> : sets time (in seconds) that the
     *    program pauses after PetscDrawPause() has been called
     *    (0 is default, -1 implies until user input).
     */
    void printVecASCII( 
            const Vec &b,           ///< [in] Petsc Vec 1xn vector to print to screen
            char const name[]="printVecASCII.txt", ///< [in] filename to write to
            const PetscViewerFormat format=PETSC_VIEWER_DEFAULT ///< [in] format for viewer
            );
}

#endif
