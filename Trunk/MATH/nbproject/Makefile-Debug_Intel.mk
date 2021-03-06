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
CND_CONF=Debug_Intel
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/LinearAlgebra.o \
	${OBJECTDIR}/src/MC_Iterative_Jacobi.o \
	${OBJECTDIR}/src/ParallelLinearAlgebra.o \
	${OBJECTDIR}/src/Point2D.o \
	${OBJECTDIR}/src/Point3D.o \
	${OBJECTDIR}/src/Vector2D.o \
	${OBJECTDIR}/src/Vector3D.o


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
LDLIBSOPTIONS=-Wl,-rpath,../UTILS/dist/Debug_Intel/Intel-Linux-x86 -L../UTILS/dist/Debug_Intel/Intel-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}: ../UTILS/dist/Debug_Intel/Intel-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	icpc -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/LinearAlgebra.o: src/LinearAlgebra.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -I.. -I../UTILS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/LinearAlgebra.o src/LinearAlgebra.c

${OBJECTDIR}/src/MC_Iterative_Jacobi.o: src/MC_Iterative_Jacobi.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -I.. -I../UTILS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MC_Iterative_Jacobi.o src/MC_Iterative_Jacobi.c

${OBJECTDIR}/src/ParallelLinearAlgebra.o: src/ParallelLinearAlgebra.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ParallelLinearAlgebra.o src/ParallelLinearAlgebra.cpp

${OBJECTDIR}/src/Point2D.o: src/Point2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Point2D.o src/Point2D.cpp

${OBJECTDIR}/src/Point3D.o: src/Point3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Point3D.o src/Point3D.cpp

${OBJECTDIR}/src/Vector2D.o: src/Vector2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vector2D.o src/Vector2D.cpp

${OBJECTDIR}/src/Vector3D.o: src/Vector3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vector3D.o src/Vector3D.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_Intel

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_Intel clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
