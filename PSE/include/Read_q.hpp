#ifndef READ_Q_H
#define READ_Q_H
#include <petscksp.h>

namespace PSE
{
    /** 
     * \brief Read q from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            Vec &evec,    ///< [in,out] eigen vector to return after reading file (already initialized)
            PetscScalar y[], ///< [out] y values to return after reading
            const int ny,    ///< [in] length of y
            PetscScalar z[], ///< [out] z values to return after reading
            const int nz,                  ///< [in] size of vector z vector
            PetscScalar &alpha, ///< [out] eigen value alpha to read
            char const buff[]="../OrrSommerfeld_and_primitive/uvwP_201"///< [in] filename of binary to read (output of python script
            );
    /** 
     * \brief Read q from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            Vec &output,    ///< [in,out] vector to return after reading file (already initialized)
            const int n,                  ///< [in] int size of vector
            char const buff[]="./src/tofile.dat"///< [in] filename of binary to read (output of python script
            );
    /** 
     * \brief Read y from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            PetscScalar output[],    ///< [out] vector to return after reading file
            const int n,                  ///< [in] int size of vector
            char const buff[]="./src/tofile.dat"///< [in] filename of binary to read (output of python script
            );
    /** 
     * \brief Read complex scalar from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            PetscScalar &output,    ///< [out] complex scalar to return after reading file
            char const buff[]="./src/tofile.dat"///< [in] filename of binary to read (output of python script
            );
}
#endif
