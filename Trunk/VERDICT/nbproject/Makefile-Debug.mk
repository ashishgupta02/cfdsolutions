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
	${OBJECTDIR}/src/V_QuadMetric.o \
	${OBJECTDIR}/src/V_TetMetric.o \
	${OBJECTDIR}/src/V_PyramidMetric.o \
	${OBJECTDIR}/src/V_GaussIntegration.o \
	${OBJECTDIR}/src/VerdictVector.o \
	${OBJECTDIR}/src/V_EdgeMetric.o \
	${OBJECTDIR}/src/V_HexMetric.o \
	${OBJECTDIR}/src/V_WedgeMetric.o \
	${OBJECTDIR}/src/V_KnifeMetric.o \
	${OBJECTDIR}/src/V_TriMetric.o

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
	${MAKE}  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/libVERDICT.so

dist/Debug/GNU-Linux-x86/libVERDICT.so: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.cc} -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVERDICT.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/V_QuadMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_QuadMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_QuadMetric.o src/V_QuadMetric.cpp

${OBJECTDIR}/src/V_TetMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_TetMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_TetMetric.o src/V_TetMetric.cpp

${OBJECTDIR}/src/V_PyramidMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_PyramidMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_PyramidMetric.o src/V_PyramidMetric.cpp

${OBJECTDIR}/src/V_GaussIntegration.o: nbproject/Makefile-${CND_CONF}.mk src/V_GaussIntegration.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_GaussIntegration.o src/V_GaussIntegration.cpp

${OBJECTDIR}/src/VerdictVector.o: nbproject/Makefile-${CND_CONF}.mk src/VerdictVector.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/VerdictVector.o src/VerdictVector.cpp

${OBJECTDIR}/src/V_EdgeMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_EdgeMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_EdgeMetric.o src/V_EdgeMetric.cpp

${OBJECTDIR}/src/V_HexMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_HexMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_HexMetric.o src/V_HexMetric.cpp

${OBJECTDIR}/src/V_WedgeMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_WedgeMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_WedgeMetric.o src/V_WedgeMetric.cpp

${OBJECTDIR}/src/V_KnifeMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_KnifeMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_KnifeMetric.o src/V_KnifeMetric.cpp

${OBJECTDIR}/src/V_TriMetric.o: nbproject/Makefile-${CND_CONF}.mk src/V_TriMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_TriMetric.o src/V_TriMetric.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/libVERDICT.so

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
