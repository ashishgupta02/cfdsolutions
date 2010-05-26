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
	${OBJECTDIR}/src/Command_Line.o \
	${OBJECTDIR}/src/Interface_MPI_Lib.o \
	${OBJECTDIR}/src/Interface_MPI.o \
	${OBJECTDIR}/src/Check_Scanf.o \
	${OBJECTDIR}/src/Linked_List.o \
	${OBJECTDIR}/src/Interface_MPI_Set.o \
	${OBJECTDIR}/src/Check_String.o \
	${OBJECTDIR}/src/Utils.o \
	${OBJECTDIR}/src/List.o \
	${OBJECTDIR}/src/Check_Malloc.o \
	${OBJECTDIR}/src/Formatting_Control.o \
	${OBJECTDIR}/src/Warn_Error.o \
	${OBJECTDIR}/src/MemUtils.o \
	${OBJECTDIR}/src/MUtils.o \
	${OBJECTDIR}/src/Logging.o \
	${OBJECTDIR}/src/Machine_Parameters.o \
	${OBJECTDIR}/src/Log.o \
	${OBJECTDIR}/src/Stopwatch.o \
	${OBJECTDIR}/src/Timing.o \
	${OBJECTDIR}/src/MathError.o

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
	${MAKE}  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/libUTILS.so

dist/Debug/GNU-Linux-x86/libUTILS.so: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.cc} -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libUTILS.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/Command_Line.o: nbproject/Makefile-${CND_CONF}.mk src/Command_Line.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Command_Line.o src/Command_Line.c

${OBJECTDIR}/src/Interface_MPI_Lib.o: nbproject/Makefile-${CND_CONF}.mk src/Interface_MPI_Lib.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Interface_MPI_Lib.o src/Interface_MPI_Lib.c

${OBJECTDIR}/src/Interface_MPI.o: nbproject/Makefile-${CND_CONF}.mk src/Interface_MPI.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Interface_MPI.o src/Interface_MPI.c

${OBJECTDIR}/src/Check_Scanf.o: nbproject/Makefile-${CND_CONF}.mk src/Check_Scanf.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Check_Scanf.o src/Check_Scanf.c

${OBJECTDIR}/src/Linked_List.o: nbproject/Makefile-${CND_CONF}.mk src/Linked_List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Linked_List.o src/Linked_List.cpp

${OBJECTDIR}/src/Interface_MPI_Set.o: nbproject/Makefile-${CND_CONF}.mk src/Interface_MPI_Set.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Interface_MPI_Set.o src/Interface_MPI_Set.c

${OBJECTDIR}/src/Check_String.o: nbproject/Makefile-${CND_CONF}.mk src/Check_String.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Check_String.o src/Check_String.c

${OBJECTDIR}/src/Utils.o: nbproject/Makefile-${CND_CONF}.mk src/Utils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Utils.o src/Utils.c

${OBJECTDIR}/src/List.o: nbproject/Makefile-${CND_CONF}.mk src/List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/List.o src/List.cpp

${OBJECTDIR}/src/Check_Malloc.o: nbproject/Makefile-${CND_CONF}.mk src/Check_Malloc.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Check_Malloc.o src/Check_Malloc.c

${OBJECTDIR}/src/Formatting_Control.o: nbproject/Makefile-${CND_CONF}.mk src/Formatting_Control.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Formatting_Control.o src/Formatting_Control.c

${OBJECTDIR}/src/Warn_Error.o: nbproject/Makefile-${CND_CONF}.mk src/Warn_Error.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Warn_Error.o src/Warn_Error.c

${OBJECTDIR}/src/MemUtils.o: nbproject/Makefile-${CND_CONF}.mk src/MemUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MemUtils.o src/MemUtils.cpp

${OBJECTDIR}/src/MUtils.o: nbproject/Makefile-${CND_CONF}.mk src/MUtils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MUtils.o src/MUtils.c

${OBJECTDIR}/src/Logging.o: nbproject/Makefile-${CND_CONF}.mk src/Logging.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Logging.o src/Logging.c

${OBJECTDIR}/src/Machine_Parameters.o: nbproject/Makefile-${CND_CONF}.mk src/Machine_Parameters.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Machine_Parameters.o src/Machine_Parameters.c

${OBJECTDIR}/src/Log.o: nbproject/Makefile-${CND_CONF}.mk src/Log.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Log.o src/Log.c

${OBJECTDIR}/src/Stopwatch.o: nbproject/Makefile-${CND_CONF}.mk src/Stopwatch.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Stopwatch.o src/Stopwatch.c

${OBJECTDIR}/src/Timing.o: nbproject/Makefile-${CND_CONF}.mk src/Timing.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Timing.o src/Timing.c

${OBJECTDIR}/src/MathError.o: nbproject/Makefile-${CND_CONF}.mk src/MathError.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MathError.o src/MathError.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/libUTILS.so

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
