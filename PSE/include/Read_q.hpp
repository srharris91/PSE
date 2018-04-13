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
            Vec &output,    ///< vector to return after reading file
            int n,                  ///< int size of vector
            char const buff[]="./src/tofile.dat"///< filename of binary to read (output of python script
            );
}
#endif