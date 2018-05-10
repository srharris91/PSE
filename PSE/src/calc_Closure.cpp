#include "calc_Closure.hpp"
#include "Init_Vec.hpp"
#include "set_Vec.hpp"
#include "trapz.hpp"
namespace PSE
{
    int calc_Closure(
            const Vec &q,
            Vec &qp1,
            const PetscInt &ny,
            const PetscInt &nz,
            PetscScalar &I,
            const PetscScalar &Deltay,
            const PetscScalar &Deltaz
            ){ 
        // size of q and qp1
        PetscInt n=4*ny*nz;
        // calc dq = q_{i+1} - q_i = qp1-q
        Vec dq;
        Init_Vec(dq,n);
        VecCopy(qp1,dq);    // dq=qp1
        VecAXPY(dq,-1.,q);  // dq=dq-q_i
        set_Vec(dq);        // assemble dq

        // take conjugate of qp1
        Vec qpconj;
        Init_Vec(qpconj,n);
        VecCopy(qp1,qpconj);
        VecConjugate(qpconj);
        set_Vec(qpconj);    // assemble qpconj

        // pointwise  multiply qpconj and dq
        Vec qtqx;
        Init_Vec(qtqx,n);
        VecPointwiseMult(qtqx,qpconj,dq);
        set_Vec(qtqx);

        // integrate
        trapz(qtqx,ny,nz,I,Deltay,Deltaz);

        // free memory
        PetscErrorCode ierr;
        ierr = VecDestroy(&dq); CHKERRQ(ierr);
        ierr = VecDestroy(&qpconj); CHKERRQ(ierr);
        ierr = VecDestroy(&qtqx); CHKERRQ(ierr);


        return 0;
    }
}
