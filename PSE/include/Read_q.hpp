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
            Vec &evec,    ///< eigen vector to return after reading file
            PetscScalar y[], ///< y values to return after reading
            const int ny,    ///< length of y
            PetscScalar z[], ///< z values to return after reading
            const int nz,                  ///< size of vector z vector
            PetscScalar &alpha, ///< eigen value alpha to read
            char const buff[]="../OrrSommerfeld_and_primitive/uvwP_201"///< filename of binary to read (output of python script
            );
    /** 
     * \brief Read q from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            Vec &output,    ///< vector to return after reading file
            const int n,                  ///< int size of vector
            char const buff[]="./src/tofile.dat"///< filename of binary to read (output of python script
            );
    /** 
     * \brief Read y from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            PetscScalar output[],    ///< vector to return after reading file
            const int n,                  ///< int size of vector
            char const buff[]="./src/tofile.dat"///< filename of binary to read (output of python script
            );
    /** 
     * \brief Read complex scalar from a tofile output from python
     * \return 0 if successful
     */
    int Read_q(
            PetscScalar &output,    ///< complex scalar to return after reading file
            char const buff[]="./src/tofile.dat"///< filename of binary to read (output of python script
            );
}
#endif
