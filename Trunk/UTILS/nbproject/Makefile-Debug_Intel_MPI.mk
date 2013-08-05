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
CC=icc
CCC=icpc
CXX=icpc
FC=ifort
AS=as

# Macros
CND_PLATFORM=Intel-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug_Intel_MPI
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Check_Malloc.o \
	${OBJECTDIR}/src/Check_Scanf.o \
	${OBJECTDIR}/src/Check_String.o \
	${OBJECTDIR}/src/Command_Line.o \
	${OBJECTDIR}/src/DoubleList.o \
	${OBJECTDIR}/src/Formatting_Control.o \
	${OBJECTDIR}/src/Interface_MPI.o \
	${OBJECTDIR}/src/Interface_MPI_Lib.o \
	${OBJECTDIR}/src/Interface_MPI_Set.o \
	${OBJECTDIR}/src/Linked_List.o \
	${OBJECTDIR}/src/List.o \
	${OBJECTDIR}/src/Log.o \
	${OBJECTDIR}/src/Logging.o \
	${OBJECTDIR}/src/MUtils.o \
	${OBJECTDIR}/src/Machine_Parameters.o \
	${OBJECTDIR}/src/MathError.o \
	${OBJECTDIR}/src/MemUtils.o \
	${OBJECTDIR}/src/Stopwatch.o \
	${OBJECTDIR}/src/Timing.o \
	${OBJECTDIR}/src/Trim_Utils.o \
	${OBJECTDIR}/src/Utils.o \
	${OBJECTDIR}/src/Warn_Error.o


# C Compiler Flags
CFLAGS=-fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000

# CC Compiler Flags
CCFLAGS=-fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000
CXXFLAGS=-fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000

# Fortran Compiler Flags
FFLAGS=-fast

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/apps/HPC_LIBS/lib

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libUTILS.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libUTILS.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	icpc -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libUTILS.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -lmpi -shared -fPIC

${OBJECTDIR}/src/Check_Malloc.o: src/Check_Malloc.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Check_Malloc.o src/Check_Malloc.c

${OBJECTDIR}/src/Check_Scanf.o: src/Check_Scanf.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Check_Scanf.o src/Check_Scanf.c

${OBJECTDIR}/src/Check_String.o: src/Check_String.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Check_String.o src/Check_String.c

${OBJECTDIR}/src/Command_Line.o: src/Command_Line.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Command_Line.o src/Command_Line.c

${OBJECTDIR}/src/DoubleList.o: src/DoubleList.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DoubleList.o src/DoubleList.cpp

${OBJECTDIR}/src/Formatting_Control.o: src/Formatting_Control.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Formatting_Control.o src/Formatting_Control.c

${OBJECTDIR}/src/Interface_MPI.o: src/Interface_MPI.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Interface_MPI.o src/Interface_MPI.c

${OBJECTDIR}/src/Interface_MPI_Lib.o: src/Interface_MPI_Lib.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Interface_MPI_Lib.o src/Interface_MPI_Lib.c

${OBJECTDIR}/src/Interface_MPI_Set.o: src/Interface_MPI_Set.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Interface_MPI_Set.o src/Interface_MPI_Set.c

${OBJECTDIR}/src/Linked_List.o: src/Linked_List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Linked_List.o src/Linked_List.cpp

${OBJECTDIR}/src/List.o: src/List.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/List.o src/List.cpp

${OBJECTDIR}/src/Log.o: src/Log.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Log.o src/Log.c

${OBJECTDIR}/src/Logging.o: src/Logging.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Logging.o src/Logging.c

${OBJECTDIR}/src/MUtils.o: src/MUtils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MUtils.o src/MUtils.c

${OBJECTDIR}/src/Machine_Parameters.o: src/Machine_Parameters.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Machine_Parameters.o src/Machine_Parameters.c

${OBJECTDIR}/src/MathError.o: src/MathError.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MathError.o src/MathError.c

${OBJECTDIR}/src/MemUtils.o: src/MemUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MemUtils.o src/MemUtils.cpp

${OBJECTDIR}/src/Stopwatch.o: src/Stopwatch.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Stopwatch.o src/Stopwatch.c

${OBJECTDIR}/src/Timing.o: src/Timing.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Timing.o src/Timing.c

${OBJECTDIR}/src/Trim_Utils.o: src/Trim_Utils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Trim_Utils.o src/Trim_Utils.c

${OBJECTDIR}/src/Utils.o: src/Utils.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Utils.o src/Utils.c

${OBJECTDIR}/src/Warn_Error.o: src/Warn_Error.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DHAVE_MPI -I.. -Iinclude -I/apps/HPC_LIBS/include -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Warn_Error.o src/Warn_Error.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libUTILS.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
