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
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/GEOMETRY/src/Vector2D.o \
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
	${OBJECTDIR}/GEOMETRY/src/Point2D.o \
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
LDLIBSOPTIONS=-lm -Wl,-rpath ../UTILS/dist/Release/GNU-Linux-x86 -L../UTILS/dist/Release/GNU-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d: ../UTILS/dist/Release/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -lstdc++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/GEOMETRY/src/Vector2D.o: GEOMETRY/src/Vector2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Vector2D.o GEOMETRY/src/Vector2D.cpp

${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o: MESHIO/src/TriangleMeshDB.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o MESHIO/src/TriangleMeshDB.cpp

${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o: MESHER/src/TriangleMeshSLK.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMeshSLK.o MESHER/src/TriangleMeshSLK.cpp

${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o: ADAPTREFINE/src/Edge2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/ADAPTREFINE/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ADAPTREFINE/src/Edge2D.o ADAPTREFINE/src/Edge2D.cpp

${OBJECTDIR}/MESHER/src/Bowyer_Watson.o: MESHER/src/Bowyer_Watson.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Bowyer_Watson.o MESHER/src/Bowyer_Watson.cpp

${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o: WINSLOW/src/TriangleWinslowSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOW/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o WINSLOW/src/TriangleWinslowSmoother.cpp

${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o: WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOWVV/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp

${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o: LESMOOTH/src/TriangleLinearElasticSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/LESMOOTH/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o LESMOOTH/src/TriangleLinearElasticSmoother.cpp

${OBJECTDIR}/MESHER/src/TriangleMesh.o: MESHER/src/TriangleMesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesh.o MESHER/src/TriangleMesh.cpp

${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o: ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp 
	${MKDIR} -p ${OBJECTDIR}/ADAPTREFINE/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp

${OBJECTDIR}/MESHER/src/TriangleMesher2D.o: MESHER/src/TriangleMesher2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesher2D.o MESHER/src/TriangleMesher2D.cpp

${OBJECTDIR}/MESHER/src/Grid_Utils.o: MESHER/src/Grid_Utils.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Grid_Utils.o MESHER/src/Grid_Utils.cpp

${OBJECTDIR}/GEOMETRY/src/Point2D.o: GEOMETRY/src/Point2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Point2D.o GEOMETRY/src/Point2D.cpp

${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o: OPTIMESH/src/TriangleMeshOptimizer.cpp 
	${MKDIR} -p ${OBJECTDIR}/OPTIMESH/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -IMESHER/include -IOPTIMESH/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -IGEOMETRY/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o OPTIMESH/src/TriangleMeshOptimizer.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
