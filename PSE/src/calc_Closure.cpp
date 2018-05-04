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
            const PetscScalar &Deltax
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

        // integrate
        I=0;
        // for each variable integrate in y
        // now integrate each in z
        
        // sub part
        Vec qsub;
        Init_Vec(qsub,ny);
        PetscScalar trap_value;
        //for (int zi=0; zi<nz; ++zi){
        int zi=0; // currently integrate over just one zi plane.... TODO
        for (int vari=0; vari<4; ++vari){
            PetscInt row_eq = (4*zi + vari)*ny;
            set_Vec(qtqx,row_eq,row_eq+ny,qsub);
            trapz(qsub,ny,trap_value,Deltax);
            I += trap_value;
        }
        //}

        // free memory
        PetscErrorCode ierr;
        ierr = VecDestroy(&dq); CHKERRQ(ierr);
        ierr = VecDestroy(&qpconj); CHKERRQ(ierr);
        ierr = VecDestroy(&qtqx); CHKERRQ(ierr);
        ierr = VecDestroy(&qsub); CHKERRQ(ierr);


        return 0;
    }
}
