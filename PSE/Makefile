PETSC_DIR = Submodules/petsc
INCLUDE_DIR      = ./include
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         = -std=c++14 -I${INCLUDE_DIR}
FPPFLAGS         =
NP               = 3
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
	${MPIEXEC} -n ${NP} ./PSE -pc_type none -ksp_gmres_cgs_refinement_type refine_always -draw_pause -1
	-@#-on_error_attach_debugger
	-@#${MPIEXEC} -n ${NP} ./PSE -print_x6 -pc_type none -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always #> PSE_1.tmp 2>&1;
	-@# ${RM} -f PSE_1.tmp
