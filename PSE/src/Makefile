PETSC_DIR = ../Submodules/petsc
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
NP               = 1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex11: ex11.o 
	-${CLINKER} -o ex11 ex11.o  ${PETSC_KSP_LIB}
	${RM} ex11.o

ex1: ex1.o 
	-${CLINKER} -o ex1 ex1.o  ${PETSC_KSP_LIB}
	${RM} ex1.o
ex2: ex2.o 
	-${CLINKER} -o ex2 ex2.o  ${PETSC_KSP_LIB}
	${RM} ex2.o
#----------------------------------------------------------------------------
runex11:
	-@${MPIEXEC} -n 1 ./ex11 -n 6 -norandom -pc_type none -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex11_1.tmp 2>&1; \
	   #${RM} -f ex11_1.tmp
runex1:
		-@${MPIEXEC} -n 1 ./ex1 -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex1_1.tmp 2>&1;	  \
		  #${RM} -f ex1_1.tmp

runex2:
		-@${MPIEXEC} -n 1 ./ex2 -ksp_monitor_short -m 5 -n 5 -ksp_gmres_cgs_refinement_type refine_always > ex2_1.tmp 2>&1; \
		  #${RM} -f ex2_1.tmp

