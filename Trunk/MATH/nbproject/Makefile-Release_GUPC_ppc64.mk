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
CND_PLATFORM=GUPC_ppc64-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release_GUPC_ppc64
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
CFLAGS=-mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000

# CC Compiler Flags
CCFLAGS=-mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000
CXXFLAGS=-mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000

# Fortran Compiler Flags
FFLAGS=-mcpu=powerpc64 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fdefault-double-8 -fdefault-real-8 -I./include

# Assembler Flags
ASFLAGS=-march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include

# Link Libraries and Options
LDLIBSOPTIONS=-Wl,-rpath,../UTILS/dist/Release_GUPC_ppc64/GUPC_ppc64-Linux-x86 -L../UTILS/dist/Release_GUPC_ppc64/GUPC_ppc64-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}: ../UTILS/dist/Release_GUPC_ppc64/GUPC_ppc64-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -mcpu=powerpc64 -shared -fPIC

${OBJECTDIR}/src/LinearAlgebra.o: src/LinearAlgebra.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -I.. -I../UTILS/include -Iinclude -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/LinearAlgebra.o src/LinearAlgebra.c

${OBJECTDIR}/src/MC_Iterative_Jacobi.o: src/MC_Iterative_Jacobi.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -I.. -I../UTILS/include -Iinclude -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MC_Iterative_Jacobi.o src/MC_Iterative_Jacobi.c

${OBJECTDIR}/src/ParallelLinearAlgebra.o: src/ParallelLinearAlgebra.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -I../UTILS/include -Iinclude -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ParallelLinearAlgebra.o src/ParallelLinearAlgebra.cpp

${OBJECTDIR}/src/Point2D.o: src/Point2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -I../UTILS/include -Iinclude -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Point2D.o src/Point2D.cpp

${OBJECTDIR}/src/Point3D.o: src/Point3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -I../UTILS/include -Iinclude -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Point3D.o src/Point3D.cpp

${OBJECTDIR}/src/Vector2D.o: src/Vector2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -I../UTILS/include -Iinclude -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vector2D.o src/Vector2D.cpp

${OBJECTDIR}/src/Vector3D.o: src/Vector3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -I../UTILS/include -Iinclude -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vector3D.o src/Vector3D.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release_GUPC_ppc64

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release_GUPC_ppc64 clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
