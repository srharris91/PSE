#include "calc_Closure.hpp"
#include "Init_Vec.hpp"
namespace PSE
{
    int calc_Closure(
            const Vec &q,
            Vec &qp1,
            const PetscInt &ny,
            PetscScalar &I,
            const PetscScalar &Deltax=2
            ){
        Vec dq;
        Init_Vec(dq,ny);




        return 0;
    }
}
