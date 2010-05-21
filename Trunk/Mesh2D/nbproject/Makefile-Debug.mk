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
CND_CONF=Debug
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/GEOMETRY/src/Vector2D.o \
	${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o \
	${OBJECTDIR}/MESHER/src/Bowyer_Watson.o \
	${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o \
	${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o \
	${OBJECTDIR}/GEOMETRY/src/Point2D.o \
	${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o \
	${OBJECTDIR}/MESHER/src/TriangleMesher2D.o \
	${OBJECTDIR}/MESHER/src/Grid_Utils.o \
	${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o \
	${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o \
	${OBJECTDIR}/MESHER/src/TriangleMesh.o \
	${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o

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
LDLIBSOPTIONS=-lm -Wl,-rpath ../UTILS/dist/Debug/GNU-Linux-x86 -L../UTILS/dist/Debug/GNU-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/mesh2d

dist/Debug/GNU-Linux-x86/mesh2d: ../UTILS/dist/Debug/GNU-Linux-x86/libUTILS.so

dist/Debug/GNU-Linux-x86/mesh2d: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	g++ -lstdc++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/GEOMETRY/src/Vector2D.o: nbproject/Makefile-${CND_CONF}.mk GEOMETRY/src/Vector2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Vector2D.o GEOMETRY/src/Vector2D.cpp

${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o: nbproject/Makefile-${CND_CONF}.mk ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp 
	${MKDIR} -p ${OBJECTDIR}/ADAPTREFINE/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp

${OBJECTDIR}/MESHER/src/Bowyer_Watson.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/Bowyer_Watson.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Bowyer_Watson.o MESHER/src/Bowyer_Watson.cpp

${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/TriangleMeshSLK.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o MESHER/src/TriangleMeshSLK.cpp

${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o: nbproject/Makefile-${CND_CONF}.mk WINSLOW/src/TriangleWinslowSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOW/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o WINSLOW/src/TriangleWinslowSmoother.cpp

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o: nbproject/Makefile-${CND_CONF}.mk ADAPTREFINE/src/Edge2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/ADAPTREFINE/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o ADAPTREFINE/src/Edge2D.cpp

${OBJECTDIR}/GEOMETRY/src/Point2D.o: nbproject/Makefile-${CND_CONF}.mk GEOMETRY/src/Point2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Point2D.o GEOMETRY/src/Point2D.cpp

${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/TriangleMeshDB.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o MESHIO/src/TriangleMeshDB.cpp

${OBJECTDIR}/MESHER/src/TriangleMesher2D.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/TriangleMesher2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesher2D.o MESHER/src/TriangleMesher2D.cpp

${OBJECTDIR}/MESHER/src/Grid_Utils.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/Grid_Utils.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Grid_Utils.o MESHER/src/Grid_Utils.cpp

${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o: nbproject/Makefile-${CND_CONF}.mk OPTIMESH/src/TriangleMeshOptimizer.cpp 
	${MKDIR} -p ${OBJECTDIR}/OPTIMESH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o OPTIMESH/src/TriangleMeshOptimizer.cpp

${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o: nbproject/Makefile-${CND_CONF}.mk LESMOOTH/src/TriangleLinearElasticSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/LESMOOTH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o LESMOOTH/src/TriangleLinearElasticSmoother.cpp

${OBJECTDIR}/MESHER/src/TriangleMesh.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/TriangleMesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesh.o MESHER/src/TriangleMesh.cpp

${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o: nbproject/Makefile-${CND_CONF}.mk WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOWVV/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/mesh2d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
