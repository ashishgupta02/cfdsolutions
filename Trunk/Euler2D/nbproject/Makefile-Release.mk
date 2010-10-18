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
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/DESIGN/src/Euler2D_Design.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o \
	${OBJECTDIR}/MESH/src/Euler2D_Mesh.o \
	${OBJECTDIR}/MATCOMP/src/MC_Iterative_Jacobi.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o \
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
LDLIBSOPTIONS=-lm -Wl,-rpath ../UTILS/dist/Release/GNU-Linux-x86 -L../UTILS/dist/Release/GNU-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/euler2d

dist/Release/GNU-Linux-x86/euler2d: ../UTILS/dist/Release/GNU-Linux-x86/libUTILS.so

dist/Release/GNU-Linux-x86/euler2d: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o: SOLVER/src/Euler2D_Solver_VanLeer.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o SOLVER/src/Euler2D_Solver_VanLeer.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/DESIGN/src/Euler2D_Design.o: DESIGN/src/Euler2D_Design.cpp 
	${MKDIR} -p ${OBJECTDIR}/DESIGN/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/DESIGN/src/Euler2D_Design.o DESIGN/src/Euler2D_Design.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o: SOLVER/src/Euler2D_Solver_StegerWarming.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o SOLVER/src/Euler2D_Solver_StegerWarming.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o: SOLVER/src/Euler2D_Solver_Roe.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o SOLVER/src/Euler2D_Solver_Roe.cpp

${OBJECTDIR}/MESH/src/Euler2D_Mesh.o: MESH/src/Euler2D_Mesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESH/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESH/src/Euler2D_Mesh.o MESH/src/Euler2D_Mesh.cpp

${OBJECTDIR}/MATCOMP/src/MC_Iterative_Jacobi.o: MATCOMP/src/MC_Iterative_Jacobi.c 
	${MKDIR} -p ${OBJECTDIR}/MATCOMP/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MATCOMP/src/MC_Iterative_Jacobi.o MATCOMP/src/MC_Iterative_Jacobi.c

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o: SOLVER/src/Euler2D_Solver_Osher.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o SOLVER/src/Euler2D_Solver_Osher.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o: SOLVER/src/Euler2D_Solver_AUSM.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o SOLVER/src/Euler2D_Solver_AUSM.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_LDFSS.o: SOLVER/src/Euler2D_Solver_LDFSS.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMATCOMP/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_LDFSS.o SOLVER/src/Euler2D_Solver_LDFSS.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/euler2d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
