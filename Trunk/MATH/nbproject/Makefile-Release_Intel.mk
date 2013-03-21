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
CCC=icc
CXX=icc
FC=ifort
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Release_Intel
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Point3D.o \
	${OBJECTDIR}/src/Vector2D.o \
	${OBJECTDIR}/src/MC_Iterative_Jacobi.o \
	${OBJECTDIR}/src/ParallelLinearAlgebra.o \
	${OBJECTDIR}/src/Vector3D.o \
	${OBJECTDIR}/src/LinearAlgebra.o \
	${OBJECTDIR}/src/Point2D.o


# C Compiler Flags
CFLAGS=-axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer

# CC Compiler Flags
CCFLAGS=-axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer
CXXFLAGS=-axSSE4.2 -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/home/agupta/00-Applications/Lakesite/HPC_LIBS/lib -Wl,-rpath ../UTILS/dist/Release_Intel/GNU-Linux-x86 -L../UTILS/dist/Release_Intel/GNU-Linux-x86 -lUTILS

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.so: ../UTILS/dist/Release_Intel/GNU-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.so: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	icc -lmpi -shared -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.so -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/Point3D.o: src/Point3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DHAVE_MPI -I.. -I../UTILS/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Point3D.o src/Point3D.cpp

${OBJECTDIR}/src/Vector2D.o: src/Vector2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DHAVE_MPI -I.. -I../UTILS/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vector2D.o src/Vector2D.cpp

${OBJECTDIR}/src/MC_Iterative_Jacobi.o: src/MC_Iterative_Jacobi.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -DHAVE_MPI -I.. -I../UTILS/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MC_Iterative_Jacobi.o src/MC_Iterative_Jacobi.c

${OBJECTDIR}/src/ParallelLinearAlgebra.o: src/ParallelLinearAlgebra.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DHAVE_MPI -I.. -I../UTILS/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ParallelLinearAlgebra.o src/ParallelLinearAlgebra.cpp

${OBJECTDIR}/src/Vector3D.o: src/Vector3D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DHAVE_MPI -I.. -I../UTILS/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vector3D.o src/Vector3D.cpp

${OBJECTDIR}/src/LinearAlgebra.o: src/LinearAlgebra.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.c) -O2 -Wall -DHAVE_MPI -I.. -I../UTILS/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/LinearAlgebra.o src/LinearAlgebra.c

${OBJECTDIR}/src/Point2D.o: src/Point2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -Wall -DHAVE_MPI -I.. -I../UTILS/include -Iinclude -I/home/agupta/00-Applications/Lakesite/HPC_LIBS/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Point2D.o src/Point2D.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release_Intel

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMATH.so

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Release_Intel clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
