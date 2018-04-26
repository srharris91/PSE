#include "set_Euler_Advance.hpp"

namespace PSE{
    PetscInt set_Euler_Advance(
            const PetscScalar &hx,
            Mat &A,
            Mat &B
            ){
        // set A and b for output
        MatAXPY(A,1./hx,B,DIFFERENT_NONZERO_PATTERN); // A+= 1./hx * B
        //set_Mat(A);
        MatScale(B,1./hx);  // B*=1./hx
        //MatMult(B,qn,b);    // b=B*qn


        return 0;
    }

}
