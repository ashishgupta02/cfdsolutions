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
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o \
	${OBJECTDIR}/MESH/src/Euler2D_Mesh.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o \
	${OBJECTDIR}/MATCOMP/src/MC_Iterative_Jacobi.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o \
	${OBJECTDIR}/UTILS/src/MUtils.o \
	${OBJECTDIR}/UTILS/src/Utils.o \
	${OBJECTDIR}/UTILS/src/List.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_LDFSS.o

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
LDLIBSOPTIONS=-lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/euler2d

dist/Release/GNU-Linux-x86/euler2d: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o: nbproject/Makefile-${CND_CONF}.mk SOLVER/src/Euler2D_Solver_Roe.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o SOLVER/src/Euler2D_Solver_Roe.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o: nbproject/Makefile-${CND_CONF}.mk SOLVER/src/Euler2D_Solver_Osher.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o SOLVER/src/Euler2D_Solver_Osher.cpp

${OBJECTDIR}/MESH/src/Euler2D_Mesh.o: nbproject/Makefile-${CND_CONF}.mk MESH/src/Euler2D_Mesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESH/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESH/src/Euler2D_Mesh.o MESH/src/Euler2D_Mesh.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o: nbproject/Makefile-${CND_CONF}.mk SOLVER/src/Euler2D_Solver_AUSM.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o SOLVER/src/Euler2D_Solver_AUSM.cpp

${OBJECTDIR}/MATCOMP/src/MC_Iterative_Jacobi.o: nbproject/Makefile-${CND_CONF}.mk MATCOMP/src/MC_Iterative_Jacobi.c 
	${MKDIR} -p ${OBJECTDIR}/MATCOMP/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MATCOMP/src/MC_Iterative_Jacobi.o MATCOMP/src/MC_Iterative_Jacobi.c

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o: nbproject/Makefile-${CND_CONF}.mk SOLVER/src/Euler2D_Solver_StegerWarming.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o SOLVER/src/Euler2D_Solver_StegerWarming.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o: nbproject/Makefile-${CND_CONF}.mk SOLVER/src/Euler2D_Solver_VanLeer.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o SOLVER/src/Euler2D_Solver_VanLeer.cpp

${OBJECTDIR}/UTILS/src/MUtils.o: nbproject/Makefile-${CND_CONF}.mk UTILS/src/MUtils.c 
	${MKDIR} -p ${OBJECTDIR}/UTILS/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/UTILS/src/MUtils.o UTILS/src/MUtils.c

${OBJECTDIR}/UTILS/src/Utils.o: nbproject/Makefile-${CND_CONF}.mk UTILS/src/Utils.c 
	${MKDIR} -p ${OBJECTDIR}/UTILS/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/UTILS/src/Utils.o UTILS/src/Utils.c

${OBJECTDIR}/UTILS/src/List.o: nbproject/Makefile-${CND_CONF}.mk UTILS/src/List.cpp 
	${MKDIR} -p ${OBJECTDIR}/UTILS/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/UTILS/src/List.o UTILS/src/List.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_LDFSS.o: nbproject/Makefile-${CND_CONF}.mk SOLVER/src/Euler2D_Solver_LDFSS.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -IUTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_LDFSS.o SOLVER/src/Euler2D_Solver_LDFSS.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/euler2d

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc