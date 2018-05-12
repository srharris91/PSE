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
            PSE();      ///< constructor with no arguments (set default values for everything)
            void init();     ///< initialize matrices and vectors
            ~PSE();     ///< destructor
            PetscInt destroy();  ///< free memory

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
            Vec b,      ///< RHS vector for Ax=b PSE step
                q,      ///< \f$ q_i \f$
                qp1;    ///< \f$ q_{i+1} \f$
            PetscScalar Re,     ///< Reynolds number
                        rho,    ///< density of fluid
                        alpha,  ///< \f$ \alpha \f$ spatial eigenvalue (x-direction)
                        m,      ///< \f$m\f$
                        omega,  ///< \f$ \omega \f$ temporal eigenvalue
                        hx,     ///< distance to march in x-direction
                        xsteps; ///< number of steps to take
            PetscScalar *y,     ///< y vector (each processor has whole vector)
                        *z;     ///< z vector (each processor has whole vector)
            PetscInt ny,        ///< number of grid points in the y-direction
                     nz,        ///< number of grid points in the z-direction
                     dim,       ///< dimension of large matrices \f$ ny*nz \f$
                     order;     ///< order of accuracy
            // flags
            PetscBool reduce_wall_order=PETSC_TRUE, ///< do we want to reduce the order of the derivative at the wall? (only in y-direction because z is periodic)
                      flag_set_D=PETSC_FALSE,     ///< flag to set derivative operators
                      flag_set_I=PETSC_FALSE,     ///< flag to set identity matrix
                      flag_base_flow=PETSC_FALSE, ///< flag to set base flow U, Uy, etc.
                      flag_init=PETSC_FALSE;      ///< flag to initialize all matrices and vectors on MPI
            PetscErrorCode ierr;            ///< flag to catch error from Petsc functions

    };
}
#endif
