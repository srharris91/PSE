#ifndef PSE_H
#define PSE_H
/**@file
 * \mainpage PSE Solver
 * \author Shaun Harris (<A HREF="https://srharris91.github.io/" TARGET="_top">https://srharris91.github.io/</A>)\n
 * Copyright (C) 2018\n
 * \section Info Information
 * This code will compute the Parabolic Stability Equation (PSE) for various flows.  
 * This will read in the Orr-Sommerfeld solution from a python script, and will march it forward in \f$x\f$.
 * \image html Diagram_closure.png "Solution Procedure" width=10cm
 *
 * All of the code is contained in the namespace PSE.  
 *
 * The PSE equation is defined here
 * \f[
 * A \hat{q}_{n,m} + 
 * B \frac{\partial \hat{q}_{n,m}}{\partial y} + 
 * C \frac{\partial^2 \hat{q}_{n,m}}{\partial y^2} + 
 * D \frac{\partial \hat{q}_{n,m}}{\partial x} = 
 * F
 *
 * \f]
 * \section compileinfo Petsc Compile Info
 * Make sure to compile petsc with options
 * \code{.sh}
 * ./configure --with-scalar-type=complex --with-precision=double
 *  \endcode
 *
 *
 * \section DownloadSourceCode Download
 * Download PSE Solver source code using following GitHub link:\n
 * <A HREF="https://github.com/srharris91/PSE.git" TARGET="_top">https://github.com/srharris91/PSE.git</A>
 *
 * \section PETSc PETSc documentation
 * This documentation contains the doxygen generated files for the PETSc library.  However, PETSc was written and used with a different documentation generator known as Sowing. http://wgropp.cs.illinois.edu/projects/software/sowing/
 *
 * The full documentation for PETSc can be found here: <A HREF="http://www.mcs.anl.gov/petsc/petsc-current/docs/manual.pdf" TARGET="_top">http://www.mcs.anl.gov/petsc/petsc-current/docs/manual.pdf</A>
 * and here: <A HREF="https://www.mcs.anl.gov/petsc/index.html" TARGET="_top">https://www.mcs.anl.gov/petsc/index.html</A>
 *
 * */
#include <petscksp.h>
#include <petscsys.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "base_flow.hpp"
#include "Ax_b.hpp"
#include "Read_q.hpp"
#include "print.hpp"
#include "get_D_Coeffs.hpp"
#include "Init_Vec.hpp"
#include "Init_Mat.hpp"
#include "set_D.hpp"
#include "set_Mat.hpp"
#include "set_Vec.hpp"
#include "set_A_and_B_zi.hpp"
#include "set_A_and_B.hpp"
#include "set_BCs.hpp"
#include "set_b.hpp"
#include "set_Euler_Advance.hpp"
#include "trapz.hpp"
#include "calc_Closure.hpp"
#include "update_Closure.hpp"

#endif
