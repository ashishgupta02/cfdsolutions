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
	${OBJECTDIR}/src/MAIN.o \
	${OBJECTDIR}/src/MESHQAZone.o \
	${OBJECTDIR}/src/MESHQACFDDB.o \
	${OBJECTDIR}/src/MESHQACellInfo.o \
	${OBJECTDIR}/src/MESHQARoot.o \
	${OBJECTDIR}/src/MESHQAEnums.o \
	${OBJECTDIR}/src/MESHQAElemset.o \
	${OBJECTDIR}/src/MESHQABase.o \
	${OBJECTDIR}/src/MESHQACellObject.o \
	${OBJECTDIR}/src/MESHQACellQualityVector.o \
	${OBJECTDIR}/src/MESHQAQualityMetrics.o

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
LDLIBSOPTIONS=-Wl,-rpath ../VERDICT/dist/Debug/GNU-Linux-x86 -L../VERDICT/dist/Debug/GNU-Linux-x86 -lVERDICT

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/libMESHQA.so

dist/Debug/GNU-Linux-x86/libMESHQA.so: ../VERDICT/dist/Debug/GNU-Linux-x86/libVERDICT.so

dist/Debug/GNU-Linux-x86/libMESHQA.so: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.cc} -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMESHQA.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/MAIN.o: nbproject/Makefile-${CND_CONF}.mk src/MAIN.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MAIN.o src/MAIN.cpp

${OBJECTDIR}/src/MESHQAZone.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQAZone.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQAZone.o src/MESHQAZone.cpp

${OBJECTDIR}/src/MESHQACFDDB.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQACFDDB.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQACFDDB.o src/MESHQACFDDB.cpp

${OBJECTDIR}/src/MESHQACellInfo.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQACellInfo.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQACellInfo.o src/MESHQACellInfo.cpp

${OBJECTDIR}/src/MESHQARoot.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQARoot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQARoot.o src/MESHQARoot.cpp

${OBJECTDIR}/src/MESHQAEnums.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQAEnums.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQAEnums.o src/MESHQAEnums.cpp

${OBJECTDIR}/src/MESHQAElemset.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQAElemset.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQAElemset.o src/MESHQAElemset.cpp

${OBJECTDIR}/src/MESHQABase.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQABase.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQABase.o src/MESHQABase.cpp

${OBJECTDIR}/src/MESHQACellObject.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQACellObject.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQACellObject.o src/MESHQACellObject.cpp

${OBJECTDIR}/src/MESHQACellQualityVector.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQACellQualityVector.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQACellQualityVector.o src/MESHQACellQualityVector.cpp

${OBJECTDIR}/src/MESHQAQualityMetrics.o: nbproject/Makefile-${CND_CONF}.mk src/MESHQAQualityMetrics.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DLINUX_AMD64 -DMESHQA_SHARED_LIB -Iinclude -I../UTILS/include -I../CFDDB/include -I../VERDICT/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MESHQAQualityMetrics.o src/MESHQAQualityMetrics.cpp

# Subprojects
.build-subprojects:
	cd ../VERDICT && ${MAKE}  -f Makefile CONF=Debug
	cd ../VERDICT && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/libMESHQA.so

# Subprojects
.clean-subprojects:
	cd ../VERDICT && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../VERDICT && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
