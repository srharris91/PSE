#ifndef PRINT_H
#define PRINT_H

#include <petscksp.h>
namespace PSE
{
    /** Print PetscScalar variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printScalar(
            PetscScalar x[],    ///< PetscScalar array to print to screen
            int n=1,              ///< size of scalar array to print
            char const name[]="x"     ///< name of variable to output default to 'x'
            ); 

    /** Print Vec from PETSc type variable to screen
     * \return ierr from PetscErrorCode of PetscPrintf
     */
    PetscErrorCode printVec(
            Vec &x,    ///< Vec array to print to screen
            int n=1,              ///< size of scalar array to print
            char const name[]="x"     ///< name of variable to output default to 'x'
            ); 
}

#endif
