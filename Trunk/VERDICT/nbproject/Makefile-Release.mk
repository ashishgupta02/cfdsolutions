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
	${OBJECTDIR}/src/V_WedgeMetric.o \
	${OBJECTDIR}/src/V_PyramidMetric.o \
	${OBJECTDIR}/src/VerdictVector.o \
	${OBJECTDIR}/src/V_TriMetric.o \
	${OBJECTDIR}/src/V_QuadMetric.o \
	${OBJECTDIR}/src/V_GaussIntegration.o \
	${OBJECTDIR}/src/V_EdgeMetric.o \
	${OBJECTDIR}/src/V_HexMetric.o \
	${OBJECTDIR}/src/V_KnifeMetric.o \
	${OBJECTDIR}/src/V_TetMetric.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVERDICT.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVERDICT.so: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVERDICT.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/V_WedgeMetric.o: src/V_WedgeMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_WedgeMetric.o src/V_WedgeMetric.cpp

${OBJECTDIR}/src/V_PyramidMetric.o: src/V_PyramidMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_PyramidMetric.o src/V_PyramidMetric.cpp

${OBJECTDIR}/src/VerdictVector.o: src/VerdictVector.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/VerdictVector.o src/VerdictVector.cpp

${OBJECTDIR}/src/V_TriMetric.o: src/V_TriMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_TriMetric.o src/V_TriMetric.cpp

${OBJECTDIR}/src/V_QuadMetric.o: src/V_QuadMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_QuadMetric.o src/V_QuadMetric.cpp

${OBJECTDIR}/src/V_GaussIntegration.o: src/V_GaussIntegration.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_GaussIntegration.o src/V_GaussIntegration.cpp

${OBJECTDIR}/src/V_EdgeMetric.o: src/V_EdgeMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_EdgeMetric.o src/V_EdgeMetric.cpp

${OBJECTDIR}/src/V_HexMetric.o: src/V_HexMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_HexMetric.o src/V_HexMetric.cpp

${OBJECTDIR}/src/V_KnifeMetric.o: src/V_KnifeMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_KnifeMetric.o src/V_KnifeMetric.cpp

${OBJECTDIR}/src/V_TetMetric.o: src/V_TetMetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DBUILD_SHARED_LIBS -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/V_TetMetric.o src/V_TetMetric.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVERDICT.so

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
