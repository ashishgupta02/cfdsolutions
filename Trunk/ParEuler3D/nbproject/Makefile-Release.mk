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
CC=mpicc
CCC=mpic++
CXX=mpic++
FC=mpif90
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
	${OBJECTDIR}/MESHIO/src/Mesh3D_IO.o \
	${OBJECTDIR}/MESHIO/src/Node3D.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/GEOMETRY/src/Vector2D.o \
	${OBJECTDIR}/MESHIO/src/Cell3D.o \
	${OBJECTDIR}/CommMPI/src/CommMPI.o \
	${OBJECTDIR}/MESHIO/src/BC3D.o \
	${OBJECTDIR}/GEOMETRY/src/Vector3D.o \
	${OBJECTDIR}/MESHIO/src/Ghost3D.o \
	${OBJECTDIR}/MESHIO/src/Grid3D.o \
	${OBJECTDIR}/GEOMETRY/src/Point2D.o \
	${OBJECTDIR}/GEOMETRY/src/Point3D.o \
	${OBJECTDIR}/MESHIO/src/Face3D.o \
	${OBJECTDIR}/Commons.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer
CXXFLAGS=-Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L../../../../00-Applications/HPC_LIBS/lib -Wl,-rpath ../UTILS/dist/Release/GNU-Linux-x86 -L../UTILS/dist/Release/GNU-Linux-x86 -lUTILS -Wl,-rpath ../CFDDB/dist/Release/GNU-Linux-x86 -L../CFDDB/dist/Release/GNU-Linux-x86 -lCFDDB -lcgns -lnetcdf -lparmetis -lmetis

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pareuler3d

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pareuler3d: ../UTILS/dist/Release/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pareuler3d: ../CFDDB/dist/Release/GNU-Linux-x86/libCFDDB.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pareuler3d: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	mpic++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pareuler3d ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/MESHIO/src/Mesh3D_IO.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/Mesh3D_IO.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/Mesh3D_IO.o MESHIO/src/Mesh3D_IO.cpp

${OBJECTDIR}/MESHIO/src/Node3D.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/Node3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/Node3D.o MESHIO/src/Node3D.cpp

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/GEOMETRY/src/Vector2D.o: nbproject/Makefile-${CND_CONF}.mk GEOMETRY/src/Vector2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Vector2D.o GEOMETRY/src/Vector2D.cpp

${OBJECTDIR}/MESHIO/src/Cell3D.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/Cell3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/Cell3D.o MESHIO/src/Cell3D.cpp

${OBJECTDIR}/CommMPI/src/CommMPI.o: nbproject/Makefile-${CND_CONF}.mk CommMPI/src/CommMPI.cpp 
	${MKDIR} -p ${OBJECTDIR}/CommMPI/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/CommMPI/src/CommMPI.o CommMPI/src/CommMPI.cpp

${OBJECTDIR}/MESHIO/src/BC3D.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/BC3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/BC3D.o MESHIO/src/BC3D.cpp

${OBJECTDIR}/GEOMETRY/src/Vector3D.o: nbproject/Makefile-${CND_CONF}.mk GEOMETRY/src/Vector3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Vector3D.o GEOMETRY/src/Vector3D.cpp

${OBJECTDIR}/MESHIO/src/Ghost3D.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/Ghost3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/Ghost3D.o MESHIO/src/Ghost3D.cpp

${OBJECTDIR}/MESHIO/src/Grid3D.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/Grid3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/Grid3D.o MESHIO/src/Grid3D.cpp

${OBJECTDIR}/GEOMETRY/src/Point2D.o: nbproject/Makefile-${CND_CONF}.mk GEOMETRY/src/Point2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Point2D.o GEOMETRY/src/Point2D.cpp

${OBJECTDIR}/GEOMETRY/src/Point3D.o: nbproject/Makefile-${CND_CONF}.mk GEOMETRY/src/Point3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/GEOMETRY/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GEOMETRY/src/Point3D.o GEOMETRY/src/Point3D.cpp

${OBJECTDIR}/MESHIO/src/Face3D.o: nbproject/Makefile-${CND_CONF}.mk MESHIO/src/Face3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/MESHIO/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MESHIO/src/Face3D.o MESHIO/src/Face3D.cpp

${OBJECTDIR}/Commons.o: nbproject/Makefile-${CND_CONF}.mk Commons.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -I../UTILS/include -I../CFDDB/include -I. -ICommMPI/include -IGEOMETRY/include -IMESHIO/include -ISOLVER/include -I../../../../00-Applications/HPC_LIBS/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/Commons.o Commons.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release
	cd ../CFDDB && ${MAKE}  -f Makefile CONF=Release
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release
	cd ../CFDDB && ${MAKE}  -f Makefile CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pareuler3d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean
	cd ../CFDDB && ${MAKE}  -f Makefile CONF=Release clean
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release clean
	cd ../CFDDB && ${MAKE}  -f Makefile CONF=Release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
