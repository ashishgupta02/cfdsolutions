#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=icc
CCC=icc
CXX=icc
FC=ifort
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug_Intel
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/AUSM_Fluxes.o \
	${OBJECTDIR}/src/Area_Volume.o \
	${OBJECTDIR}/src/BC.o \
	${OBJECTDIR}/src/Commons.o \
	${OBJECTDIR}/src/CompressibleUtils.o \
	${OBJECTDIR}/src/Connectivity_Maps.o \
	${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o \
	${OBJECTDIR}/src/DebugSolver.o \
	${OBJECTDIR}/src/Gradient.o \
	${OBJECTDIR}/src/HLLC_Fluxes.o \
	${OBJECTDIR}/src/HigherOrderReconstructQ.o \
	${OBJECTDIR}/src/JST_Fluxes.o \
	${OBJECTDIR}/src/Jacobian.o \
	${OBJECTDIR}/src/LDFSS_Fluxes.o \
	${OBJECTDIR}/src/Limiters.o \
	${OBJECTDIR}/src/Material.o \
	${OBJECTDIR}/src/MeshIO.o \
	${OBJECTDIR}/src/Osher_Fluxes.o \
	${OBJECTDIR}/src/RMSIO.o \
	${OBJECTDIR}/src/Residual.o \
	${OBJECTDIR}/src/Residual_Smoothing.o \
	${OBJECTDIR}/src/RestartIO.o \
	${OBJECTDIR}/src/Roe_EntropyFix.o \
	${OBJECTDIR}/src/Roe_Fluxes.o \
	${OBJECTDIR}/src/Roe_Jacobian.o \
	${OBJECTDIR}/src/Solver.o \
	${OBJECTDIR}/src/SolverParameters.o \
	${OBJECTDIR}/src/Solver_Steady_Explicit.o \
	${OBJECTDIR}/src/Solver_Steady_Implicit.o \
	${OBJECTDIR}/src/Solver_Unsteady_Explicit.o \
	${OBJECTDIR}/src/Solver_Unsteady_Implicit.o \
	${OBJECTDIR}/src/StegerWarming_Fluxes.o \
	${OBJECTDIR}/src/Test.o \
	${OBJECTDIR}/src/Time_Step.o \
	${OBJECTDIR}/src/VanLeer_Fluxes.o


# C Compiler Flags
CFLAGS=-axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer

# CC Compiler Flags
CCFLAGS=-axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer
CXXFLAGS=-axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/home/agupta/00-Applications/Lakesite/HPC_LIBS/lib -Wl,-rpath,../MATH/dist/Debug_Intel/GNU-Linux-x86 -L../MATH/dist/Debug_Intel/GNU-Linux-x86 -lMATH -Wl,-rpath,../UTILS/dist/Debug_Intel/GNU-Linux-x86 -L../UTILS/dist/Debug_Intel/GNU-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ../MATH/dist/Debug_Intel/GNU-Linux-x86/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ../UTILS/dist/Debug_Intel/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	icc -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d ${OBJECTFILES} ${LDLIBSOPTIONS} -lstdc++ -lmpi

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/AUSM_Fluxes.o: src/AUSM_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/AUSM_Fluxes.o src/AUSM_Fluxes.cpp

${OBJECTDIR}/src/Area_Volume.o: src/Area_Volume.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Area_Volume.o src/Area_Volume.cpp

${OBJECTDIR}/src/BC.o: src/BC.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/BC.o src/BC.cpp

${OBJECTDIR}/src/Commons.o: src/Commons.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Commons.o src/Commons.cpp

${OBJECTDIR}/src/CompressibleUtils.o: src/CompressibleUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/CompressibleUtils.o src/CompressibleUtils.cpp

${OBJECTDIR}/src/Connectivity_Maps.o: src/Connectivity_Maps.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Connectivity_Maps.o src/Connectivity_Maps.cpp

${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o: src/Cuthill_Mckee_Reorder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o src/Cuthill_Mckee_Reorder.cpp

${OBJECTDIR}/src/DebugSolver.o: src/DebugSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DebugSolver.o src/DebugSolver.cpp

${OBJECTDIR}/src/Gradient.o: src/Gradient.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Gradient.o src/Gradient.cpp

${OBJECTDIR}/src/HLLC_Fluxes.o: src/HLLC_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/HLLC_Fluxes.o src/HLLC_Fluxes.cpp

${OBJECTDIR}/src/HigherOrderReconstructQ.o: src/HigherOrderReconstructQ.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/HigherOrderReconstructQ.o src/HigherOrderReconstructQ.cpp

${OBJECTDIR}/src/JST_Fluxes.o: src/JST_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/JST_Fluxes.o src/JST_Fluxes.cpp

${OBJECTDIR}/src/Jacobian.o: src/Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Jacobian.o src/Jacobian.cpp

${OBJECTDIR}/src/LDFSS_Fluxes.o: src/LDFSS_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/LDFSS_Fluxes.o src/LDFSS_Fluxes.cpp

${OBJECTDIR}/src/Limiters.o: src/Limiters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Limiters.o src/Limiters.cpp

${OBJECTDIR}/src/Material.o: src/Material.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Material.o src/Material.cpp

${OBJECTDIR}/src/MeshIO.o: src/MeshIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MeshIO.o src/MeshIO.cpp

${OBJECTDIR}/src/Osher_Fluxes.o: src/Osher_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Osher_Fluxes.o src/Osher_Fluxes.cpp

${OBJECTDIR}/src/RMSIO.o: src/RMSIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/RMSIO.o src/RMSIO.cpp

${OBJECTDIR}/src/Residual.o: src/Residual.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Residual.o src/Residual.cpp

${OBJECTDIR}/src/Residual_Smoothing.o: src/Residual_Smoothing.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Residual_Smoothing.o src/Residual_Smoothing.cpp

${OBJECTDIR}/src/RestartIO.o: src/RestartIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/RestartIO.o src/RestartIO.cpp

${OBJECTDIR}/src/Roe_EntropyFix.o: src/Roe_EntropyFix.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_EntropyFix.o src/Roe_EntropyFix.cpp

${OBJECTDIR}/src/Roe_Fluxes.o: src/Roe_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_Fluxes.o src/Roe_Fluxes.cpp

${OBJECTDIR}/src/Roe_Jacobian.o: src/Roe_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_Jacobian.o src/Roe_Jacobian.cpp

${OBJECTDIR}/src/Solver.o: src/Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver.o src/Solver.cpp

${OBJECTDIR}/src/SolverParameters.o: src/SolverParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/SolverParameters.o src/SolverParameters.cpp

${OBJECTDIR}/src/Solver_Steady_Explicit.o: src/Solver_Steady_Explicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Steady_Explicit.o src/Solver_Steady_Explicit.cpp

${OBJECTDIR}/src/Solver_Steady_Implicit.o: src/Solver_Steady_Implicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Steady_Implicit.o src/Solver_Steady_Implicit.cpp

${OBJECTDIR}/src/Solver_Unsteady_Explicit.o: src/Solver_Unsteady_Explicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Unsteady_Explicit.o src/Solver_Unsteady_Explicit.cpp

${OBJECTDIR}/src/Solver_Unsteady_Implicit.o: src/Solver_Unsteady_Implicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Unsteady_Implicit.o src/Solver_Unsteady_Implicit.cpp

${OBJECTDIR}/src/StegerWarming_Fluxes.o: src/StegerWarming_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/StegerWarming_Fluxes.o src/StegerWarming_Fluxes.cpp

${OBJECTDIR}/src/Test.o: src/Test.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Test.o src/Test.cpp

${OBJECTDIR}/src/Time_Step.o: src/Time_Step.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Time_Step.o src/Time_Step.cpp

${OBJECTDIR}/src/VanLeer_Fluxes.o: src/VanLeer_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -I../UTILS/include -I../MATH/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/VanLeer_Fluxes.o src/VanLeer_Fluxes.cpp

# Subprojects
.build-subprojects:
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug_Intel
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_Intel

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d

# Subprojects
.clean-subprojects:
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug_Intel clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_Intel clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
