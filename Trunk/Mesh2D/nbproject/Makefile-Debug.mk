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
AS=

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
	${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o \
	${OBJECTDIR}/MESHER/src/Bowyer_Watson.o \
	${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o \
	${OBJECTDIR}/MESHER/src/TriangleMesher2D.o \
	${OBJECTDIR}/MESHER/src/Grid_Utils.o \
	${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o \
	${OBJECTDIR}/UTILS/src/MUtils.o \
	${OBJECTDIR}/UTILS/src/Utils.o \
	${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o \
	${OBJECTDIR}/UTILS/src/List.o \
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
LDLIBSOPTIONS=-L../../00-Applications/HPC_LIBS/lib -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/mesh2d

dist/Debug/GNU-Linux-x86/mesh2d: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mesh2d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o: nbproject/Makefile-${CND_CONF}.mk ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp 
	${MKDIR} -p ${OBJECTDIR}/ADAPTREFINE/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ADAPTREFINE/src/TriangleAdaptiveRefinement.o ADAPTREFINE/src/TriangleAdaptiveRefinement.cpp

${OBJECTDIR}/MESHER/src/Bowyer_Watson.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/Bowyer_Watson.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Bowyer_Watson.o MESHER/src/Bowyer_Watson.cpp

${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o: nbproject/Makefile-${CND_CONF}.mk WINSLOW/src/TriangleWinslowSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOW/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOW/src/TriangleWinslowSmoother.o WINSLOW/src/TriangleWinslowSmoother.cpp

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/TriangleMeshDB.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/TriangleMeshDB.o MESHIO/src/TriangleMeshDB.cpp

${OBJECTDIR}/MESHER/src/TriangleMesher2D.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/TriangleMesher2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesher2D.o MESHER/src/TriangleMesher2D.cpp

${OBJECTDIR}/MESHER/src/Grid_Utils.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/Grid_Utils.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/Grid_Utils.o MESHER/src/Grid_Utils.cpp

${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o: nbproject/Makefile-${CND_CONF}.mk OPTIMESH/src/TriangleMeshOptimizer.cpp 
	${MKDIR} -p ${OBJECTDIR}/OPTIMESH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/OPTIMESH/src/TriangleMeshOptimizer.o OPTIMESH/src/TriangleMeshOptimizer.cpp

${OBJECTDIR}/UTILS/src/MUtils.o: nbproject/Makefile-${CND_CONF}.mk UTILS/src/MUtils.c 
	${MKDIR} -p ${OBJECTDIR}/UTILS/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/UTILS/src/MUtils.o UTILS/src/MUtils.c

${OBJECTDIR}/UTILS/src/Utils.o: nbproject/Makefile-${CND_CONF}.mk UTILS/src/Utils.c 
	${MKDIR} -p ${OBJECTDIR}/UTILS/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/UTILS/src/Utils.o UTILS/src/Utils.c

${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o: nbproject/Makefile-${CND_CONF}.mk LESMOOTH/src/TriangleLinearElasticSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/LESMOOTH/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/LESMOOTH/src/TriangleLinearElasticSmoother.o LESMOOTH/src/TriangleLinearElasticSmoother.cpp

${OBJECTDIR}/UTILS/src/List.o: nbproject/Makefile-${CND_CONF}.mk UTILS/src/List.cpp 
	${MKDIR} -p ${OBJECTDIR}/UTILS/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/UTILS/src/List.o UTILS/src/List.cpp

${OBJECTDIR}/MESHER/src/TriangleMesh.o: nbproject/Makefile-${CND_CONF}.mk MESHER/src/TriangleMesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHER/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHER/src/TriangleMesh.o MESHER/src/TriangleMesh.cpp

${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o: nbproject/Makefile-${CND_CONF}.mk WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp 
	${MKDIR} -p ${OBJECTDIR}/WINSLOWVV/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DBUILD_LIB -IMESHER/include -IOPTIMESH/include -IUTILS/include -IWINSLOW/include -ILESMOOTH/include -IWINSLOWVV/include -IADAPTREFINE/include -IMESHIO/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.o WINSLOWVV/src/TriangleWinslowVirtualVolumeSmoother.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/mesh2d

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
