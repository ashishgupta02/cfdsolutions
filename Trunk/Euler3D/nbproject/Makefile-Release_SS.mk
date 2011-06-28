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
CC=cc
CCC=CC
CXX=CC
FC=f95
AS=as

# Macros
CND_PLATFORM=OracleSolarisStudio-Linux-x86
CND_CONF=Release_SS
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Point3D.o \
	${OBJECTDIR}/src/Limiters.o \
	${OBJECTDIR}/_ext/979082451/HigherOrderReconstructQ.o \
	${OBJECTDIR}/src/Commons.o \
	${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/RestartIO.o \
	${OBJECTDIR}/src/MC_Iterative_Jacobi.o \
	${OBJECTDIR}/src/Vector3D.o \
	${OBJECTDIR}/src/BC.o \
	${OBJECTDIR}/_ext/979082451/Roe_EntropyFix.o \
	${OBJECTDIR}/src/Time_Step.o \
	${OBJECTDIR}/src/LMRoe_Fluxes.o \
	${OBJECTDIR}/src/Gradient.o \
	${OBJECTDIR}/src/Trim_Utils.o \
	${OBJECTDIR}/src/Area_Volume.o \
	${OBJECTDIR}/src/List.o \
	${OBJECTDIR}/_ext/979082451/CompressibleUtils.o \
	${OBJECTDIR}/src/MeshIO.o \
	${OBJECTDIR}/src/Roe_Fluxes.o \
	${OBJECTDIR}/src/Solver.o \
	${OBJECTDIR}/src/Connectivity_Maps.o \
	${OBJECTDIR}/src/DebugSolver.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-Release_SS.mk dist/Release_SS/OracleSolarisStudio-Linux-x86/euler3d

dist/Release_SS/OracleSolarisStudio-Linux-x86/euler3d: ${OBJECTFILES}
	${MKDIR} -p dist/Release_SS/OracleSolarisStudio-Linux-x86
	CC -lm -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/Point3D.o: nbproject/Makefile-${CND_CONF}.mk src/Point3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Point3D.o src/Point3D.cpp

${OBJECTDIR}/src/Limiters.o: nbproject/Makefile-${CND_CONF}.mk src/Limiters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Limiters.o src/Limiters.cpp

${OBJECTDIR}/_ext/979082451/HigherOrderReconstructQ.o: nbproject/Makefile-${CND_CONF}.mk /home/agupta/NetBeansProjects/CFDSolutions/Trunk/Euler3D/src/HigherOrderReconstructQ.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/979082451
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/_ext/979082451/HigherOrderReconstructQ.o /home/agupta/NetBeansProjects/CFDSolutions/Trunk/Euler3D/src/HigherOrderReconstructQ.cpp

${OBJECTDIR}/src/Commons.o: nbproject/Makefile-${CND_CONF}.mk src/Commons.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Commons.o src/Commons.cpp

${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o: nbproject/Makefile-${CND_CONF}.mk src/Cuthill_Mckee_Reorder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o src/Cuthill_Mckee_Reorder.cpp

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/RestartIO.o: nbproject/Makefile-${CND_CONF}.mk src/RestartIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/RestartIO.o src/RestartIO.cpp

${OBJECTDIR}/src/MC_Iterative_Jacobi.o: nbproject/Makefile-${CND_CONF}.mk src/MC_Iterative_Jacobi.c 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.c) -fast -g +w -Iinclude -o ${OBJECTDIR}/src/MC_Iterative_Jacobi.o src/MC_Iterative_Jacobi.c

${OBJECTDIR}/src/Vector3D.o: nbproject/Makefile-${CND_CONF}.mk src/Vector3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Vector3D.o src/Vector3D.cpp

${OBJECTDIR}/src/BC.o: nbproject/Makefile-${CND_CONF}.mk src/BC.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/BC.o src/BC.cpp

${OBJECTDIR}/_ext/979082451/Roe_EntropyFix.o: nbproject/Makefile-${CND_CONF}.mk /home/agupta/NetBeansProjects/CFDSolutions/Trunk/Euler3D/src/Roe_EntropyFix.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/979082451
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/_ext/979082451/Roe_EntropyFix.o /home/agupta/NetBeansProjects/CFDSolutions/Trunk/Euler3D/src/Roe_EntropyFix.cpp

${OBJECTDIR}/src/Time_Step.o: nbproject/Makefile-${CND_CONF}.mk src/Time_Step.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Time_Step.o src/Time_Step.cpp

${OBJECTDIR}/src/LMRoe_Fluxes.o: nbproject/Makefile-${CND_CONF}.mk src/LMRoe_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/LMRoe_Fluxes.o src/LMRoe_Fluxes.cpp

${OBJECTDIR}/src/Gradient.o: nbproject/Makefile-${CND_CONF}.mk src/Gradient.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Gradient.o src/Gradient.cpp

${OBJECTDIR}/src/Trim_Utils.o: nbproject/Makefile-${CND_CONF}.mk src/Trim_Utils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.c) -fast -g +w -Iinclude -o ${OBJECTDIR}/src/Trim_Utils.o src/Trim_Utils.c

${OBJECTDIR}/src/Area_Volume.o: nbproject/Makefile-${CND_CONF}.mk src/Area_Volume.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Area_Volume.o src/Area_Volume.cpp

${OBJECTDIR}/src/List.o: nbproject/Makefile-${CND_CONF}.mk src/List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/List.o src/List.cpp

${OBJECTDIR}/_ext/979082451/CompressibleUtils.o: nbproject/Makefile-${CND_CONF}.mk /home/agupta/NetBeansProjects/CFDSolutions/Trunk/Euler3D/src/CompressibleUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/979082451
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/_ext/979082451/CompressibleUtils.o /home/agupta/NetBeansProjects/CFDSolutions/Trunk/Euler3D/src/CompressibleUtils.cpp

${OBJECTDIR}/src/MeshIO.o: nbproject/Makefile-${CND_CONF}.mk src/MeshIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/MeshIO.o src/MeshIO.cpp

${OBJECTDIR}/src/Roe_Fluxes.o: nbproject/Makefile-${CND_CONF}.mk src/Roe_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Roe_Fluxes.o src/Roe_Fluxes.cpp

${OBJECTDIR}/src/Solver.o: nbproject/Makefile-${CND_CONF}.mk src/Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Solver.o src/Solver.cpp

${OBJECTDIR}/src/Connectivity_Maps.o: nbproject/Makefile-${CND_CONF}.mk src/Connectivity_Maps.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/Connectivity_Maps.o src/Connectivity_Maps.cpp

${OBJECTDIR}/src/DebugSolver.o: nbproject/Makefile-${CND_CONF}.mk src/DebugSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 +w -Iinclude -o ${OBJECTDIR}/src/DebugSolver.o src/DebugSolver.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release_SS
	${RM} dist/Release_SS/OracleSolarisStudio-Linux-x86/euler3d
	${CCADMIN} -clean

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
