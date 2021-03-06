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
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/MESH/src/Euler2D_Mesh_LinearElasticSmoother.o \
	${OBJECTDIR}/DESIGN/src/Euler2D_Design.o \
	${OBJECTDIR}/DESIGN/src/Euler2D_Design_VanLeer.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o \
	${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o \
	${OBJECTDIR}/MESH/src/Euler2D_Mesh.o \
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
LDLIBSOPTIONS=-lm -Wl,-rpath ../UTILS/dist/Debug/GNU-Linux-x86 -L../UTILS/dist/Debug/GNU-Linux-x86 -lUTILS -Wl,-rpath ../MATH/dist/Debug/GNU-Linux-x86 -L../MATH/dist/Debug/GNU-Linux-x86 -lMATH

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d: ../UTILS/dist/Debug/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d: ../MATH/dist/Debug/GNU-Linux-x86/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	g++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o: SOLVER/src/Euler2D_Solver_VanLeer.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_VanLeer.o SOLVER/src/Euler2D_Solver_VanLeer.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/MESH/src/Euler2D_Mesh_LinearElasticSmoother.o: MESH/src/Euler2D_Mesh_LinearElasticSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESH/src/Euler2D_Mesh_LinearElasticSmoother.o MESH/src/Euler2D_Mesh_LinearElasticSmoother.cpp

${OBJECTDIR}/DESIGN/src/Euler2D_Design.o: DESIGN/src/Euler2D_Design.cpp 
	${MKDIR} -p ${OBJECTDIR}/DESIGN/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/DESIGN/src/Euler2D_Design.o DESIGN/src/Euler2D_Design.cpp

${OBJECTDIR}/DESIGN/src/Euler2D_Design_VanLeer.o: DESIGN/src/Euler2D_Design_VanLeer.cpp 
	${MKDIR} -p ${OBJECTDIR}/DESIGN/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/DESIGN/src/Euler2D_Design_VanLeer.o DESIGN/src/Euler2D_Design_VanLeer.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver.o: SOLVER/src/Euler2D_Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver.o SOLVER/src/Euler2D_Solver.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o: SOLVER/src/Euler2D_Solver_StegerWarming.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_StegerWarming.o SOLVER/src/Euler2D_Solver_StegerWarming.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o: SOLVER/src/Euler2D_Solver_Roe.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Roe.o SOLVER/src/Euler2D_Solver_Roe.cpp

${OBJECTDIR}/MESH/src/Euler2D_Mesh.o: MESH/src/Euler2D_Mesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESH/src/Euler2D_Mesh.o MESH/src/Euler2D_Mesh.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o: SOLVER/src/Euler2D_Solver_Osher.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_Osher.o SOLVER/src/Euler2D_Solver_Osher.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o: SOLVER/src/Euler2D_Solver_AUSM.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_AUSM.o SOLVER/src/Euler2D_Solver_AUSM.cpp

${OBJECTDIR}/SOLVER/src/Euler2D_Solver_LDFSS.o: SOLVER/src/Euler2D_Solver_LDFSS.cpp 
	${MKDIR} -p ${OBJECTDIR}/SOLVER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DVERBOSE -DDEBUG -I.. -I../UTILS/include -I../MATH/include -IMESH/include -ISOLVER/include -IDESIGN/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SOLVER/src/Euler2D_Solver_LDFSS.o SOLVER/src/Euler2D_Solver_LDFSS.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/euler2d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
