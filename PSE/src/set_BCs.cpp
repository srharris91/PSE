#include "set_BCs.hpp"

namespace PSE
{
    PetscInt set_BCs(
            Mat &A,
            Vec &b,
            const PetscInt &ny,
            const PetscInt &nz
            ){
        for (PetscInt zi=0; zi<nz; zi++){
            // u mom bottom and top
            PetscInt i=(4*zi+0)*ny;
            MatZeroRows(A,1,&i,1.,0,0);
            i=(4*zi+1)*ny - 1;
            MatZeroRows(A,1,&i,1.,0,0);
            // v-mom bottom and top
            i=(4*zi+1)*ny;
            MatZeroRows(A,1,&i,1.,0,0);
            i=(4*zi+2)*ny-1;
            MatZeroRows(A,1,&i,1.,0,0);
            // w-mom bottom and top
            i=(4*zi+2)*ny;
            MatZeroRows(A,1,&i,1.,0,0);
            i=(4*zi+3)*ny - 1;
            MatZeroRows(A,1,&i,1.,0,0);
        }


        return 0;
    }
}
