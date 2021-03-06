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
CND_CONF=Debug_GUPC_ppc64
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/EOS.o \
	${OBJECTDIR}/src/EOS_Density_Pressure.o \
	${OBJECTDIR}/src/EOS_Density_Temperature.o \
	${OBJECTDIR}/src/EOS_IdealGas.o \
	${OBJECTDIR}/src/EOS_Internal.o \
	${OBJECTDIR}/src/EOS_NIST.o \
	${OBJECTDIR}/src/EOS_Pressure_Temperature.o \
	${OBJECTDIR}/src/Test.o


# C Compiler Flags
CFLAGS=-mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000

# CC Compiler Flags
CCFLAGS=-mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000
CXXFLAGS=-mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000

# Fortran Compiler Flags
FFLAGS=-mcpu=powerpc64 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fdefault-double-8 -fdefault-real-8 -I./include

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-Wl,-rpath,../MATH/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86 -L../MATH/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86 -lMATH -Wl,-rpath,../UTILS/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86 -L../UTILS/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86 -lUTILS -Wl,-rpath,../NISTThermo/NISTThermo/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86 -L../NISTThermo/NISTThermo/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86 -lNISTThermo

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ../MATH/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ../UTILS/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ../NISTThermo/NISTThermo/dist/Debug_GUPC_ppc64/GUPC_ppc64-Linux-x86/libNISTThermo.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -mcpu=powerpc64 -lgfortran -shared -fPIC

${OBJECTDIR}/src/EOS.o: nbproject/Makefile-${CND_CONF}.mk src/EOS.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS.o src/EOS.c

${OBJECTDIR}/src/EOS_Density_Pressure.o: nbproject/Makefile-${CND_CONF}.mk src/EOS_Density_Pressure.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Density_Pressure.o src/EOS_Density_Pressure.c

${OBJECTDIR}/src/EOS_Density_Temperature.o: nbproject/Makefile-${CND_CONF}.mk src/EOS_Density_Temperature.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Density_Temperature.o src/EOS_Density_Temperature.c

${OBJECTDIR}/src/EOS_IdealGas.o: nbproject/Makefile-${CND_CONF}.mk src/EOS_IdealGas.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_IdealGas.o src/EOS_IdealGas.c

${OBJECTDIR}/src/EOS_Internal.o: nbproject/Makefile-${CND_CONF}.mk src/EOS_Internal.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Internal.o src/EOS_Internal.c

${OBJECTDIR}/src/EOS_NIST.o: nbproject/Makefile-${CND_CONF}.mk src/EOS_NIST.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_NIST.o src/EOS_NIST.c

${OBJECTDIR}/src/EOS_Pressure_Temperature.o: nbproject/Makefile-${CND_CONF}.mk src/EOS_Pressure_Temperature.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Pressure_Temperature.o src/EOS_Pressure_Temperature.c

${OBJECTDIR}/src/Test.o: nbproject/Makefile-${CND_CONF}.mk src/Test.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -mcpu=powerpc64 -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Test.o src/Test.cpp

# Subprojects
.build-subprojects:
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug_GUPC_ppc64
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_GUPC_ppc64
	cd ../NISTThermo/NISTThermo && ${MAKE}  -f Makefile CONF=Debug_GUPC_ppc64

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug_GUPC_ppc64 clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_GUPC_ppc64 clean
	cd ../NISTThermo/NISTThermo && ${MAKE}  -f Makefile CONF=Debug_GUPC_ppc64 clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
