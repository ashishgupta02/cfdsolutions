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
	${OBJECTDIR}/src/ugridIO.o \
	${OBJECTDIR}/src/DBMANAGER.o \
	${OBJECTDIR}/src/netcdfIO.o \
	${OBJECTDIR}/src/DB.o \
	${OBJECTDIR}/src/DBUGRID.o \
	${OBJECTDIR}/src/corelib.o \
	${OBJECTDIR}/src/DBNETCDF.o \
	${OBJECTDIR}/src/DBERROR.o \
	${OBJECTDIR}/src/cgnsIO.o \
	${OBJECTDIR}/src/DBCGNS.o


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
LDLIBSOPTIONS=-L../../../../00-Applications/HPC_LIBS/intel64/lib -Wl,-rpath ../UTILS/dist/Release/GNU-Linux-x86 -L../UTILS/dist/Release/GNU-Linux-x86 -lUTILS -lcgns -lnetcdf

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCFDDB.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCFDDB.so: ../UTILS/dist/Release/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCFDDB.so: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCFDDB.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/ugridIO.o: src/ugridIO.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -Iinclude -I.. -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ugridIO.o src/ugridIO.c

${OBJECTDIR}/src/DBMANAGER.o: src/DBMANAGER.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBMANAGER.o src/DBMANAGER.cpp

${OBJECTDIR}/src/netcdfIO.o: src/netcdfIO.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -Iinclude -I.. -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/netcdfIO.o src/netcdfIO.c

${OBJECTDIR}/src/DB.o: src/DB.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DB.o src/DB.cpp

${OBJECTDIR}/src/DBUGRID.o: src/DBUGRID.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBUGRID.o src/DBUGRID.cpp

${OBJECTDIR}/src/corelib.o: src/corelib.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -Iinclude -I.. -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/corelib.o src/corelib.c

${OBJECTDIR}/src/DBNETCDF.o: src/DBNETCDF.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBNETCDF.o src/DBNETCDF.cpp

${OBJECTDIR}/src/DBERROR.o: src/DBERROR.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBERROR.o src/DBERROR.cpp

${OBJECTDIR}/src/cgnsIO.o: src/cgnsIO.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -Iinclude -I.. -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/cgnsIO.o src/cgnsIO.c

${OBJECTDIR}/src/DBCGNS.o: src/DBCGNS.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I.. -Iinclude -I../UTILS/include -I../../../../00-Applications/HPC_LIBS/intel64/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DBCGNS.o src/DBCGNS.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCFDDB.so

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
