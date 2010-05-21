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
	${OBJECTDIR}/src/Utils.o \
	${OBJECTDIR}/src/MemUtils.o \
	${OBJECTDIR}/src/List.o \
	${OBJECTDIR}/src/MUtils.o \
	${OBJECTDIR}/src/Linked_List.o

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
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/libUTILS.so

dist/Release/GNU-Linux-x86/libUTILS.so: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libUTILS.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/Utils.o: nbproject/Makefile-${CND_CONF}.mk src/Utils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Utils.o src/Utils.c

${OBJECTDIR}/src/MemUtils.o: nbproject/Makefile-${CND_CONF}.mk src/MemUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MemUtils.o src/MemUtils.cpp

${OBJECTDIR}/src/List.o: nbproject/Makefile-${CND_CONF}.mk src/List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/List.o src/List.cpp

${OBJECTDIR}/src/MUtils.o: nbproject/Makefile-${CND_CONF}.mk src/MUtils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MUtils.o src/MUtils.c

${OBJECTDIR}/src/Linked_List.o: nbproject/Makefile-${CND_CONF}.mk src/Linked_List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Linked_List.o src/Linked_List.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/libUTILS.so

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
