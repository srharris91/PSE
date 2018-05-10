#ifndef CLASS_H
#define CLASS_H
#include <petscksp.h>

namespace PSE
{
    /** 
     * \brief Class to contain various variables for PSE solution
     */
    class PSE {
        public:
            PSE();      ///< constructor with no arguments
            init();     ///< initialize matrices and vectors
            ~PSE();     ///< destructor
            destroy();  ///< free memory

            // memory storage
            Mat Dy,     ///< derivative in y direction
                Dyy,    ///< second derivative in y direction
                Dz,     ///< derivative in z direction
                Dzz,    ///< second derivative in z direction
                U,      ///< base flow \f$U\f$ velocity
                Uy,     ///< base flow \f$\frac{dU}{dy}\f$ 
                A,      ///< Matrix to Solve PSE equation in linear system (all q,dqdy,dqdz terms)
                B,      ///< Matrix to Solve PSE equation in linear system (all dqdx terms)
                I;      ///< identity matrix, useful in setting up A,B
            // flags
    };
}
#endif
