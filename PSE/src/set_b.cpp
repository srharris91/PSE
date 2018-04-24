#include "set_b.hpp"

namespace PSE
{
    PetscInt set_b(
            const Mat &B,
            const Vec &qn,
            Vec &b         
            ){
        // set A and b for output
        MatMult(B,qn,b);    // b=B*qn

        return 0;
    }
}
