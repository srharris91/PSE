PETSC_DIR = Submodules/petsc
INCLUDE_DIR      = ./include
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         = -std=c++14 -I${INCLUDE_DIR}
FPPFLAGS         =
NP               = 4
SRCFILES         = $(wildcard ./src/*cpp)
OBJFILES         = $(patsubst %.cpp,%.o,${SRCFILES})
EXECUTABLE       = PSE

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

comp_run_PSE: exPSE runPSE 

docs:
	doxygen

clean:: clean-legacy
	@${RM} PSE

exPSE: ${OBJFILES}
	-@${CLINKER} -o ${EXECUTABLE} ${OBJFILES}  ${PETSC_KSP_LIB}
	-@${RM} ${OBJFILES}

runPSE:
	${MPIEXEC} -n ${NP} ./PSE -ksp_converged_reason -ksp_monitor_true_residual -pc_type jacobi -ksp_gmres_cgs_refinement_type refine_always -draw_pause -1 
	-@#-ksp_rtol 1e-8 
	-@#-ksp_type cg
	-@#-pc_type none
	-@#-on_error_attach_debugger
	-@#${MPIEXEC} -n ${NP} ./PSE -print_x6 -pc_type none -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always #> PSE_1.tmp 2>&1;
	-@# ${RM} -f PSE_1.tmp
