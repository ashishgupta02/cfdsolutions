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
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/AUSM_Fluxes.o \
	${OBJECTDIR}/src/AUSM_Jacobian.o \
	${OBJECTDIR}/src/Area_Volume.o \
	${OBJECTDIR}/src/BC.o \
	${OBJECTDIR}/src/Commons.o \
	${OBJECTDIR}/src/CompressibleUtils.o \
	${OBJECTDIR}/src/Connectivity_Maps.o \
	${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o \
	${OBJECTDIR}/src/DebugSolver.o \
	${OBJECTDIR}/src/FiniteDifference_Jacobian.o \
	${OBJECTDIR}/src/Fluxes.o \
	${OBJECTDIR}/src/Gradient.o \
	${OBJECTDIR}/src/HLLC_Fluxes.o \
	${OBJECTDIR}/src/HLLC_Jacobian.o \
	${OBJECTDIR}/src/HigherOrderReconstructQ.o \
	${OBJECTDIR}/src/JST_Fluxes.o \
	${OBJECTDIR}/src/JST_Jacobian.o \
	${OBJECTDIR}/src/Jacobian.o \
	${OBJECTDIR}/src/LDFSS_Fluxes.o \
	${OBJECTDIR}/src/LDFSS_Jacobian.o \
	${OBJECTDIR}/src/Limiters.o \
	${OBJECTDIR}/src/Material.o \
	${OBJECTDIR}/src/Material_Transformations.o \
	${OBJECTDIR}/src/MeshIO.o \
	${OBJECTDIR}/src/Osher_Fluxes.o \
	${OBJECTDIR}/src/Osher_Jacobian.o \
	${OBJECTDIR}/src/RMSIO.o \
	${OBJECTDIR}/src/Residual.o \
	${OBJECTDIR}/src/Residual_Smoothing.o \
	${OBJECTDIR}/src/RestartIO.o \
	${OBJECTDIR}/src/Roe_EntropyFix.o \
	${OBJECTDIR}/src/Roe_Fluxes.o \
	${OBJECTDIR}/src/Roe_Jacobian.o \
	${OBJECTDIR}/src/Solver.o \
	${OBJECTDIR}/src/SolverDB.o \
	${OBJECTDIR}/src/SolverParameters.o \
	${OBJECTDIR}/src/Solver_Steady_Explicit.o \
	${OBJECTDIR}/src/Solver_Steady_Implicit.o \
	${OBJECTDIR}/src/Solver_Unsteady_Explicit.o \
	${OBJECTDIR}/src/Solver_Unsteady_Implicit.o \
	${OBJECTDIR}/src/StegerWarming_Fluxes.o \
	${OBJECTDIR}/src/StegerWarming_Jacobian.o \
	${OBJECTDIR}/src/Test.o \
	${OBJECTDIR}/src/Time_Step.o \
	${OBJECTDIR}/src/VanLeer_Fluxes.o \
	${OBJECTDIR}/src/VanLeer_Jacobian.o


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
LDLIBSOPTIONS=-Wl,-rpath,../UTILS/dist/Debug_Intel/Intel-Linux-x86 -L../UTILS/dist/Debug_Intel/Intel-Linux-x86 -lUTILS -Wl,-rpath,../MATH/dist/Debug_Intel/Intel-Linux-x86 -L../MATH/dist/Debug_Intel/Intel-Linux-x86 -lMATH -Wl,-rpath,../EOS/dist/Debug_Intel/Intel-Linux-x86 -L../EOS/dist/Debug_Intel/Intel-Linux-x86 -lEOS -Wl,-rpath,../NISTThermo/NISTThermo/dist/Debug_Intel/Intel-Linux-x86 -L../NISTThermo/NISTThermo/dist/Debug_Intel/Intel-Linux-x86 -lNISTThermo

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ../UTILS/dist/Debug_Intel/Intel-Linux-x86/libUTILS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ../MATH/dist/Debug_Intel/Intel-Linux-x86/libMATH.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ../EOS/dist/Debug_Intel/Intel-Linux-x86/libEOS.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ../NISTThermo/NISTThermo/dist/Debug_Intel/Intel-Linux-x86/libNISTThermo.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	ifort -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d ${OBJECTFILES} ${LDLIBSOPTIONS} -nofor-main -lstdc++

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/AUSM_Fluxes.o: src/AUSM_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/AUSM_Fluxes.o src/AUSM_Fluxes.cpp

${OBJECTDIR}/src/AUSM_Jacobian.o: src/AUSM_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/AUSM_Jacobian.o src/AUSM_Jacobian.cpp

${OBJECTDIR}/src/Area_Volume.o: src/Area_Volume.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Area_Volume.o src/Area_Volume.cpp

${OBJECTDIR}/src/BC.o: src/BC.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/BC.o src/BC.cpp

${OBJECTDIR}/src/Commons.o: src/Commons.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Commons.o src/Commons.cpp

${OBJECTDIR}/src/CompressibleUtils.o: src/CompressibleUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/CompressibleUtils.o src/CompressibleUtils.cpp

${OBJECTDIR}/src/Connectivity_Maps.o: src/Connectivity_Maps.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Connectivity_Maps.o src/Connectivity_Maps.cpp

${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o: src/Cuthill_Mckee_Reorder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Cuthill_Mckee_Reorder.o src/Cuthill_Mckee_Reorder.cpp

${OBJECTDIR}/src/DebugSolver.o: src/DebugSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DebugSolver.o src/DebugSolver.cpp

${OBJECTDIR}/src/FiniteDifference_Jacobian.o: src/FiniteDifference_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/FiniteDifference_Jacobian.o src/FiniteDifference_Jacobian.cpp

${OBJECTDIR}/src/Fluxes.o: src/Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Fluxes.o src/Fluxes.cpp

${OBJECTDIR}/src/Gradient.o: src/Gradient.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Gradient.o src/Gradient.cpp

${OBJECTDIR}/src/HLLC_Fluxes.o: src/HLLC_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/HLLC_Fluxes.o src/HLLC_Fluxes.cpp

${OBJECTDIR}/src/HLLC_Jacobian.o: src/HLLC_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/HLLC_Jacobian.o src/HLLC_Jacobian.cpp

${OBJECTDIR}/src/HigherOrderReconstructQ.o: src/HigherOrderReconstructQ.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/HigherOrderReconstructQ.o src/HigherOrderReconstructQ.cpp

${OBJECTDIR}/src/JST_Fluxes.o: src/JST_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/JST_Fluxes.o src/JST_Fluxes.cpp

${OBJECTDIR}/src/JST_Jacobian.o: src/JST_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/JST_Jacobian.o src/JST_Jacobian.cpp

${OBJECTDIR}/src/Jacobian.o: src/Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Jacobian.o src/Jacobian.cpp

${OBJECTDIR}/src/LDFSS_Fluxes.o: src/LDFSS_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/LDFSS_Fluxes.o src/LDFSS_Fluxes.cpp

${OBJECTDIR}/src/LDFSS_Jacobian.o: src/LDFSS_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/LDFSS_Jacobian.o src/LDFSS_Jacobian.cpp

${OBJECTDIR}/src/Limiters.o: src/Limiters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Limiters.o src/Limiters.cpp

${OBJECTDIR}/src/Material.o: src/Material.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Material.o src/Material.cpp

${OBJECTDIR}/src/Material_Transformations.o: src/Material_Transformations.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Material_Transformations.o src/Material_Transformations.cpp

${OBJECTDIR}/src/MeshIO.o: src/MeshIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/MeshIO.o src/MeshIO.cpp

${OBJECTDIR}/src/Osher_Fluxes.o: src/Osher_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Osher_Fluxes.o src/Osher_Fluxes.cpp

${OBJECTDIR}/src/Osher_Jacobian.o: src/Osher_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Osher_Jacobian.o src/Osher_Jacobian.cpp

${OBJECTDIR}/src/RMSIO.o: src/RMSIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/RMSIO.o src/RMSIO.cpp

${OBJECTDIR}/src/Residual.o: src/Residual.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Residual.o src/Residual.cpp

${OBJECTDIR}/src/Residual_Smoothing.o: src/Residual_Smoothing.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Residual_Smoothing.o src/Residual_Smoothing.cpp

${OBJECTDIR}/src/RestartIO.o: src/RestartIO.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/RestartIO.o src/RestartIO.cpp

${OBJECTDIR}/src/Roe_EntropyFix.o: src/Roe_EntropyFix.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_EntropyFix.o src/Roe_EntropyFix.cpp

${OBJECTDIR}/src/Roe_Fluxes.o: src/Roe_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_Fluxes.o src/Roe_Fluxes.cpp

${OBJECTDIR}/src/Roe_Jacobian.o: src/Roe_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Roe_Jacobian.o src/Roe_Jacobian.cpp

${OBJECTDIR}/src/Solver.o: src/Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver.o src/Solver.cpp

${OBJECTDIR}/src/SolverDB.o: src/SolverDB.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/SolverDB.o src/SolverDB.cpp

${OBJECTDIR}/src/SolverParameters.o: src/SolverParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/SolverParameters.o src/SolverParameters.cpp

${OBJECTDIR}/src/Solver_Steady_Explicit.o: src/Solver_Steady_Explicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Steady_Explicit.o src/Solver_Steady_Explicit.cpp

${OBJECTDIR}/src/Solver_Steady_Implicit.o: src/Solver_Steady_Implicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Steady_Implicit.o src/Solver_Steady_Implicit.cpp

${OBJECTDIR}/src/Solver_Unsteady_Explicit.o: src/Solver_Unsteady_Explicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Unsteady_Explicit.o src/Solver_Unsteady_Explicit.cpp

${OBJECTDIR}/src/Solver_Unsteady_Implicit.o: src/Solver_Unsteady_Implicit.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Solver_Unsteady_Implicit.o src/Solver_Unsteady_Implicit.cpp

${OBJECTDIR}/src/StegerWarming_Fluxes.o: src/StegerWarming_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/StegerWarming_Fluxes.o src/StegerWarming_Fluxes.cpp

${OBJECTDIR}/src/StegerWarming_Jacobian.o: src/StegerWarming_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/StegerWarming_Jacobian.o src/StegerWarming_Jacobian.cpp

${OBJECTDIR}/src/Test.o: src/Test.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Test.o src/Test.cpp

${OBJECTDIR}/src/Time_Step.o: src/Time_Step.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Time_Step.o src/Time_Step.cpp

${OBJECTDIR}/src/VanLeer_Fluxes.o: src/VanLeer_Fluxes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/VanLeer_Fluxes.o src/VanLeer_Fluxes.cpp

${OBJECTDIR}/src/VanLeer_Jacobian.o: src/VanLeer_Jacobian.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -I.. -I../UTILS/include -I../MATH/include -I../EOS/include -Iinclude -fast -Wno-write-strings -fno-math-errno -fstrict-aliasing -fomit-frame-pointer -fpermissive -funroll-all-loops -finline-limit=150000 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/VanLeer_Jacobian.o src/VanLeer_Jacobian.cpp

# Subprojects
.build-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_Intel
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug_Intel
	cd ../EOS && ${MAKE}  -f Makefile CONF=Debug_Intel
	cd ../NISTThermo/NISTThermo && ${MAKE}  -f Makefile CONF=Debug_Intel

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpflow3d

# Subprojects
.clean-subprojects:
	cd ../UTILS && ${MAKE}  -f Makefile CONF=Debug_Intel clean
	cd ../MATH && ${MAKE}  -f Makefile CONF=Debug_Intel clean
	cd ../EOS && ${MAKE}  -f Makefile CONF=Debug_Intel clean
	cd ../NISTThermo/NISTThermo && ${MAKE}  -f Makefile CONF=Debug_Intel clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
