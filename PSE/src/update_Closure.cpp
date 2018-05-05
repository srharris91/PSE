#include "update_Closure.hpp"
#include "calc_Closure.hpp"
#include "Init_Vec.hpp"
#include "set_Vec.hpp"
#include "trapz.hpp"
#include "print.hpp"
namespace PSE
{
    int update_Closure(
            const Vec &q,
            Vec &qp1,
            const PetscInt &ny,
            const PetscInt &nz,
            const PetscScalar &dx,
            PetscScalar &Ialpha,
            PetscScalar &alpha,
            const PetscReal &tol,
            const PetscScalar &Deltax
            ){ 
        // initialize
        PetscInt dim = 4*ny*nz;
        PetscScalar Ialpha_orig = Ialpha;
        PetscScalar alpha_i = alpha; // original alpha from last step (used for trapz integration)

        // set q_physical
        Vec q_physical; // physical type q=\hat{q} exp(...)
        Init_Vec(q_physical,dim);
        VecCopy(qp1,q_physical);
        VecScale(q_physical,PetscExpScalar( PETSC_i*(Ialpha_orig + (dx*(alpha_i+alpha)/2.))));
        set_Vec(q_physical); // assemble final

        // check closure
        PetscScalar closure_value;
        calc_Closure(q,qp1,ny,nz,closure_value);
        printScalar(&closure_value,1,"Closure");
        
        Vec q2,qsub;
        Init_Vec(q2,dim);
        Init_Vec(qsub,ny);
        PetscInt iter=0;
        //while(
                //(
                 //PetscAbsReal(PetscRealPart(closure_value))        >= tol
                 //|| 
                 //PetscAbsReal(PetscImaginaryPart(closure_value))   >= tol
                //)
                //){
        for (int i=0; i<100; ++i){
            // update alpha
            VecCopy(qp1,q2);
            VecAbs(q2);
            VecPow(q2,2.);
            PetscScalar q2_int=0,trap_value;
            int zi=0; // currently integrate over just one zi plane.... TODO
            for (int vari=0; vari<4; ++vari){
                PetscInt row_eq = (4*zi + vari)*ny;
                set_Vec(q2,row_eq,row_eq+ny,qsub);
                trapz(qsub,ny,trap_value,Deltax);
                q2_int += trap_value;
            }
            PetscScalar delta_alpha = - ((PETSC_i/dx) * (closure_value)/(q2_int));
            alpha = alpha + delta_alpha;

            // update qp1
            VecCopy(q_physical,qp1); // qp1 = q_physical
            set_Vec(qp1);
            VecScale(qp1,PetscExpScalar( -1.*PETSC_i*(Ialpha_orig + (dx*(alpha_i+alpha)/2.))));
            set_Vec(qp1);

            // check closure terms
            calc_Closure(q,qp1,ny,nz,closure_value);

            // update iteration and print to screen
            iter++;
            PetscPrintf(PETSC_COMM_WORLD,"closure iteration %i \n",iter);
            printScalar(&alpha,1,"alpha");
            printScalar(&delta_alpha,1,"delta_alpha");
            printScalar(&closure_value,1,"Closure");
            char filename[100];
            sprintf(filename,"printVecqp1_%d.txt",i);
            PSE::printVecASCII(qp1  ,filename);
        }
        Ialpha=Ialpha_orig + (dx*(alpha_i+alpha)/2.);

        // free memory
        PetscErrorCode ierr;
        ierr = VecDestroy(&q_physical); CHKERRQ(ierr);
        ierr = VecDestroy(&q2); CHKERRQ(ierr);
        ierr = VecDestroy(&qsub); CHKERRQ(ierr);


        return 0;
    }
}
