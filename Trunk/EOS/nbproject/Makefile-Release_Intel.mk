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
CND_CONF=Release_Intel
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
CFLAGS=-fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000

# CC Compiler Flags
CCFLAGS=-fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000
CXXFLAGS=-fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000

# Fortran Compiler Flags
FFLAGS=-fast

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-Wl,-rpath,../UTILS/dist/Release_Intel/Intel-Linux-x86 -L../UTILS/dist/Release_Intel/Intel-Linux-x86 -lUTILS -Wl,-rpath,../MATH/dist/Release_Intel/Intel-Linux-x86 -L../MATH/dist/Release_Intel/Intel-Linux-x86 -lMATH -Wl,-rpath,../NISTThermo/NISTThermo/dist/Release_Intel/Intel-Linux-x86 -L../NISTThermo/NISTThermo/dist/Release_Intel/Intel-Linux-x86 -lNISTThermo

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ../UTILS/dist/Release_Intel/Intel-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ../MATH/dist/Release_Intel/Intel-Linux-x86/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ../NISTThermo/NISTThermo/dist/Release_Intel/Intel-Linux-x86/libNISTThermo.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/EOS.o: src/EOS.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS.o src/EOS.c

${OBJECTDIR}/src/EOS_Density_Pressure.o: src/EOS_Density_Pressure.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Density_Pressure.o src/EOS_Density_Pressure.c

${OBJECTDIR}/src/EOS_Density_Temperature.o: src/EOS_Density_Temperature.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Density_Temperature.o src/EOS_Density_Temperature.c

${OBJECTDIR}/src/EOS_IdealGas.o: src/EOS_IdealGas.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_IdealGas.o src/EOS_IdealGas.c

${OBJECTDIR}/src/EOS_Internal.o: src/EOS_Internal.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Internal.o src/EOS_Internal.c

${OBJECTDIR}/src/EOS_NIST.o: src/EOS_NIST.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_NIST.o src/EOS_NIST.c

${OBJECTDIR}/src/EOS_Pressure_Temperature.o: src/EOS_Pressure_Temperature.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EOS_Pressure_Temperature.o src/EOS_Pressure_Temperature.c

${OBJECTDIR}/src/Test.o: src/Test.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -D__linux -Iinclude -I../NISTThermo/NISTThermo/include -I../UTILS/include -I../MATH/include -I../ -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Test.o src/Test.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release_Intel
	cd ../MATH && ${MAKE}  -f Makefile CONF=Release_Intel
	cd ../NISTThermo/NISTThermo && ${MAKE}  -f Makefile CONF=Release_Intel

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libEOS.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release_Intel clean
	cd ../MATH && ${MAKE}  -f Makefile CONF=Release_Intel clean
	cd ../NISTThermo/NISTThermo && ${MAKE}  -f Makefile CONF=Release_Intel clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
