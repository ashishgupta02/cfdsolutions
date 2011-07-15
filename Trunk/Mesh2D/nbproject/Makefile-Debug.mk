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
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o \
	${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o \
	${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o \
	${OBJECTDIR}/MESHER/src/Bowyer_Watson.o \
	${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o \
	${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o \
	${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o \
	${OBJECTDIR}/MESHER/src/TriangleMesh.o \
	${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o \
	${OBJECTDIR}/MESHER/src/TriangleMesher2D.o \
	${OBJECTDIR}/MESHER/src/Grid_Utils.o \
	${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d: ../UTILS/dist/Debug/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d: ../MATH/dist/Debug/GNU-Linux-x86/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	g++ -lstdc++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o: MESHIO/src/TriangleMeshDB.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o MESHIO/src/TriangleMeshDB.cpp

${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o: MESHER/src/TriangleMeshSLK.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o MESHER/src/TriangleMeshSLK.cpp

${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o: ADAPTREFINE/src/Edge2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/ADAPTREFINE/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o ADAPTREFINE/src/Edge2D.cpp

${OBJECTDIR}/MESHER/src/Bowyer_Watson.o: MESHER/src/Bowyer_Watson.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Bowyer_Watson.o MESHER/src/Bowyer_Watson.cpp

${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o: WINSLOW/src/TriangleWinslowSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOW/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o WINSLOW/src/TriangleWinslowSmoother.cpp

${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o: WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOWVV/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp

${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o: LESMOOTH/src/TriangleLinearElasticSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/LESMOOTH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o LESMOOTH/src/TriangleLinearElasticSmoother.cpp

${OBJECTDIR}/MESHER/src/TriangleMesh.o: MESHER/src/TriangleMesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesh.o MESHER/src/TriangleMesh.cpp

${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o: ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp 
	${MKDIR} -p ${OBJECTDIR}/ADAPTREFINE/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp

${OBJECTDIR}/MESHER/src/TriangleMesher2D.o: MESHER/src/TriangleMesher2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesher2D.o MESHER/src/TriangleMesher2D.cpp

${OBJECTDIR}/MESHER/src/Grid_Utils.o: MESHER/src/Grid_Utils.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Grid_Utils.o MESHER/src/Grid_Utils.cpp

${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o: OPTIMESH/src/TriangleMeshOptimizer.cpp 
	${MKDIR} -p ${OBJECTDIR}/OPTIMESH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o OPTIMESH/src/TriangleMeshOptimizer.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
