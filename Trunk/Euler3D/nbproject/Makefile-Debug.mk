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
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Limiters.o \
	${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o \
	${OBJECTDIR}/src/Commons.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/RestartIO.o \
	${OBJECTDIR}/src/BC.o \
	${OBJECTDIR}/src/Jacobian.o \
	${OBJECTDIR}/src/Time_Step.o \
	${OBJECTDIR}/src/Roe_EntropyFix.o \
	${OBJECTDIR}/src/Gradient.o \
	${OBJECTDIR}/src/CompressibleUtils.o \
	${OBJECTDIR}/src/Area_Volume.o \
	${OBJECTDIR}/src/HigherOrderReconstructQ.o \
	${OBJECTDIR}/src/Roe_Jacobian.o \
	${OBJECTDIR}/src/MeshIO.o \
	${OBJECTDIR}/src/Roe_Fluxes.o \
	${OBJECTDIR}/src/Solver.o \
	${OBJECTDIR}/src/Connectivity_Maps.o \
	${OBJECTDIR}/src/Residual.o \
	${OBJECTDIR}/src/DebugSolver.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer
CXXFLAGS=-Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-Wl,-rpath ../UTILS/dist/Debug/GNU-Linux-x86 -L../UTILS/dist/Debug/GNU-Linux-x86 -lUTILS -Wl,-rpath ../MATH/dist/Debug/GNU-Linux-x86 -L../MATH/dist/Debug/GNU-Linux-x86 -lMATH

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d: ../UTILS/dist/Debug/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d: ../MATH/dist/Debug/GNU-Linux-x86/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	g++ -lstdc++ -lm -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/Limiters.o: src/Limiters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Limiters.o src/Limiters.cpp

${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o: src/Cuthill_Mckee_Reorder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o src/Cuthill_Mckee_Reorder.cpp

${OBJECTDIR}/src/Commons.o: src/Commons.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Commons.o src/Commons.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/RestartIO.o: src/RestartIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/RestartIO.o src/RestartIO.cpp

${OBJECTDIR}/src/BC.o: src/BC.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/BC.o src/BC.cpp

${OBJECTDIR}/src/Jacobian.o: src/Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Jacobian.o src/Jacobian.cpp

${OBJECTDIR}/src/Time_Step.o: src/Time_Step.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Time_Step.o src/Time_Step.cpp

${OBJECTDIR}/src/Roe_EntropyFix.o: src/Roe_EntropyFix.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_EntropyFix.o src/Roe_EntropyFix.cpp

${OBJECTDIR}/src/Gradient.o: src/Gradient.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Gradient.o src/Gradient.cpp

${OBJECTDIR}/src/CompressibleUtils.o: src/CompressibleUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/CompressibleUtils.o src/CompressibleUtils.cpp

${OBJECTDIR}/src/Area_Volume.o: src/Area_Volume.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Area_Volume.o src/Area_Volume.cpp

${OBJECTDIR}/src/HigherOrderReconstructQ.o: src/HigherOrderReconstructQ.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/HigherOrderReconstructQ.o src/HigherOrderReconstructQ.cpp

${OBJECTDIR}/src/Roe_Jacobian.o: src/Roe_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_Jacobian.o src/Roe_Jacobian.cpp

${OBJECTDIR}/src/MeshIO.o: src/MeshIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MeshIO.o src/MeshIO.cpp

${OBJECTDIR}/src/Roe_Fluxes.o: src/Roe_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_Fluxes.o src/Roe_Fluxes.cpp

${OBJECTDIR}/src/Solver.o: src/Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver.o src/Solver.cpp

${OBJECTDIR}/src/Connectivity_Maps.o: src/Connectivity_Maps.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Connectivity_Maps.o src/Connectivity_Maps.cpp

${OBJECTDIR}/src/Residual.o: src/Residual.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Residual.o src/Residual.cpp

${OBJECTDIR}/src/DebugSolver.o: src/DebugSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DebugSolver.o src/DebugSolver.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
