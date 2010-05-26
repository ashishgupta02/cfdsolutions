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
	${OBJECTDIR}/src/DBNETCDF.o \
	${OBJECTDIR}/src/DBERROR.o \
	${OBJECTDIR}/src/DBCGNS.o \
	${OBJECTDIR}/src/DB.o \
	${OBJECTDIR}/src/cgnsIO.o \
	${OBJECTDIR}/src/DBMANAGER.o \
	${OBJECTDIR}/src/corelib.o \
	${OBJECTDIR}/src/netcdfIO.o

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
LDLIBSOPTIONS=-L../../../../00-Applications/HPC_LIBS/lib -lcgns -lnetcdf -Wl,-rpath ../UTILS/dist/Debug/GNU-Linux-x86 -L../UTILS/dist/Debug/GNU-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/libCFDDB.so

dist/Debug/GNU-Linux-x86/libCFDDB.so: ../UTILS/dist/Debug/GNU-Linux-x86/libUTILS.so

dist/Debug/GNU-Linux-x86/libCFDDB.so: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.cc} -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCFDDB.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/DBNETCDF.o: nbproject/Makefile-${CND_CONF}.mk src/DBNETCDF.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBNETCDF.o src/DBNETCDF.cpp

${OBJECTDIR}/src/DBERROR.o: nbproject/Makefile-${CND_CONF}.mk src/DBERROR.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBERROR.o src/DBERROR.cpp

${OBJECTDIR}/src/DBCGNS.o: nbproject/Makefile-${CND_CONF}.mk src/DBCGNS.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBCGNS.o src/DBCGNS.cpp

${OBJECTDIR}/src/DB.o: nbproject/Makefile-${CND_CONF}.mk src/DB.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DB.o src/DB.cpp

${OBJECTDIR}/src/cgnsIO.o: nbproject/Makefile-${CND_CONF}.mk src/cgnsIO.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/cgnsIO.o src/cgnsIO.c

${OBJECTDIR}/src/DBMANAGER.o: nbproject/Makefile-${CND_CONF}.mk src/DBMANAGER.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBMANAGER.o src/DBMANAGER.cpp

${OBJECTDIR}/src/corelib.o: nbproject/Makefile-${CND_CONF}.mk src/corelib.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/corelib.o src/corelib.c

${OBJECTDIR}/src/netcdfIO.o: nbproject/Makefile-${CND_CONF}.mk src/netcdfIO.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -g -Wall -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/netcdfIO.o src/netcdfIO.c

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/libCFDDB.so

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc