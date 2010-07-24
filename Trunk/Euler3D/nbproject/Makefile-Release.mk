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
CND_CONF=Release
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Point3D.o \
	${OBJECTDIR}/src/Commons.o \
	${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/MC_Iterative_Jacobi.o \
	${OBJECTDIR}/src/Vector3D.o \
	${OBJECTDIR}/src/BC.o \
	${OBJECTDIR}/src/Time_Step.o \
	${OBJECTDIR}/src/Trim_Utils.o \
	${OBJECTDIR}/src/Area_Volume.o \
	${OBJECTDIR}/src/List.o \
	${OBJECTDIR}/src/MeshIO.o \
	${OBJECTDIR}/src/Roe_Fluxes.o \
	${OBJECTDIR}/src/Solver.o \
	${OBJECTDIR}/src/Connectivity_Maps.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-Wno-write-strings
CXXFLAGS=-Wno-write-strings

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/euler3d

dist/Release/GNU-Linux-x86/euler3d: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler3d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/Point3D.o: src/Point3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Point3D.o src/Point3D.cpp

${OBJECTDIR}/src/Commons.o: src/Commons.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Commons.o src/Commons.cpp

${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o: src/Cuthill_Mckee_Reorder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o src/Cuthill_Mckee_Reorder.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/MC_Iterative_Jacobi.o: src/MC_Iterative_Jacobi.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MC_Iterative_Jacobi.o src/MC_Iterative_Jacobi.c

${OBJECTDIR}/src/Vector3D.o: src/Vector3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vector3D.o src/Vector3D.cpp

${OBJECTDIR}/src/BC.o: src/BC.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/BC.o src/BC.cpp

${OBJECTDIR}/src/Time_Step.o: src/Time_Step.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Time_Step.o src/Time_Step.cpp

${OBJECTDIR}/src/Trim_Utils.o: src/Trim_Utils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Trim_Utils.o src/Trim_Utils.c

${OBJECTDIR}/src/Area_Volume.o: src/Area_Volume.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Area_Volume.o src/Area_Volume.cpp

${OBJECTDIR}/src/List.o: src/List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/List.o src/List.cpp

${OBJECTDIR}/src/MeshIO.o: src/MeshIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MeshIO.o src/MeshIO.cpp

${OBJECTDIR}/src/Roe_Fluxes.o: src/Roe_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_Fluxes.o src/Roe_Fluxes.cpp

${OBJECTDIR}/src/Solver.o: src/Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver.o src/Solver.cpp

${OBJECTDIR}/src/Connectivity_Maps.o: src/Connectivity_Maps.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Connectivity_Maps.o src/Connectivity_Maps.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/euler3d

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
