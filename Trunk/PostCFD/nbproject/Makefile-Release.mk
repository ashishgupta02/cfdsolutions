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
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Algebra.o \
	${OBJECTDIR}/src/Boundary_2D.o \
	${OBJECTDIR}/src/CGNSIO.o \
	${OBJECTDIR}/src/CGNSPatch.o \
	${OBJECTDIR}/src/CellPointOperation_2D.o \
	${OBJECTDIR}/src/ColorCode_2D.o \
	${OBJECTDIR}/src/ColorPlot_2D.o \
	${OBJECTDIR}/src/ContourPlot_2D.o \
	${OBJECTDIR}/src/Convert.o \
	${OBJECTDIR}/src/DrawMesh_2D.o \
	${OBJECTDIR}/src/Error.o \
	${OBJECTDIR}/src/FluidProperties.o \
	${OBJECTDIR}/src/Gradient_2D.o \
	${OBJECTDIR}/src/Graphics_2D.o \
	${OBJECTDIR}/src/Interpolation_2D.o \
	${OBJECTDIR}/src/Main.o \
	${OBJECTDIR}/src/PPForceMoments_2D.o \
	${OBJECTDIR}/src/PPInitializeSolution_2D.o \
	${OBJECTDIR}/src/PPOptions_2D.o \
	${OBJECTDIR}/src/PPScalarFlowVariables_2D.o \
	${OBJECTDIR}/src/PPSolutionGradients_2D.o \
	${OBJECTDIR}/src/PPTensorFlowVariables_2D.o \
	${OBJECTDIR}/src/PPVectorFlowVariables_2D.o \
	${OBJECTDIR}/src/PostAnalysis_2D.o \
	${OBJECTDIR}/src/PostProcessingInitialize_2D.o \
	${OBJECTDIR}/src/PostProcessing_2D.o \
	${OBJECTDIR}/src/PrePostProcessing_2D.o \
	${OBJECTDIR}/src/Read_CGNS.o \
	${OBJECTDIR}/src/ReferenceQuantities.o \
	${OBJECTDIR}/src/ScalarTopology_2D.o \
	${OBJECTDIR}/src/Stream_2D.o \
	${OBJECTDIR}/src/Streamline_2D.o \
	${OBJECTDIR}/src/Vector.o \
	${OBJECTDIR}/src/VectorPlot_2D.o \
	${OBJECTDIR}/src/VectorTopology_2D.o


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
LDLIBSOPTIONS=-lm -lglut -lGL -lGLU lib/sources/lib/libcgns.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/postcfd

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/postcfd: lib/sources/lib/libcgns.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/postcfd: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/postcfd ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/src/Algebra.o: src/Algebra.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Algebra.o src/Algebra.c

${OBJECTDIR}/src/Boundary_2D.o: src/Boundary_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Boundary_2D.o src/Boundary_2D.c

${OBJECTDIR}/src/CGNSIO.o: src/CGNSIO.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/CGNSIO.o src/CGNSIO.c

${OBJECTDIR}/src/CGNSPatch.o: src/CGNSPatch.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/CGNSPatch.o src/CGNSPatch.c

${OBJECTDIR}/src/CellPointOperation_2D.o: src/CellPointOperation_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/CellPointOperation_2D.o src/CellPointOperation_2D.c

${OBJECTDIR}/src/ColorCode_2D.o: src/ColorCode_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ColorCode_2D.o src/ColorCode_2D.c

${OBJECTDIR}/src/ColorPlot_2D.o: src/ColorPlot_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ColorPlot_2D.o src/ColorPlot_2D.c

${OBJECTDIR}/src/ContourPlot_2D.o: src/ContourPlot_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ContourPlot_2D.o src/ContourPlot_2D.c

${OBJECTDIR}/src/Convert.o: src/Convert.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Convert.o src/Convert.c

${OBJECTDIR}/src/DrawMesh_2D.o: src/DrawMesh_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/DrawMesh_2D.o src/DrawMesh_2D.c

${OBJECTDIR}/src/Error.o: src/Error.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Error.o src/Error.c

${OBJECTDIR}/src/FluidProperties.o: src/FluidProperties.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/FluidProperties.o src/FluidProperties.c

${OBJECTDIR}/src/Gradient_2D.o: src/Gradient_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Gradient_2D.o src/Gradient_2D.c

${OBJECTDIR}/src/Graphics_2D.o: src/Graphics_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Graphics_2D.o src/Graphics_2D.c

${OBJECTDIR}/src/Interpolation_2D.o: src/Interpolation_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Interpolation_2D.o src/Interpolation_2D.c

${OBJECTDIR}/src/Main.o: src/Main.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Main.o src/Main.c

${OBJECTDIR}/src/PPForceMoments_2D.o: src/PPForceMoments_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PPForceMoments_2D.o src/PPForceMoments_2D.c

${OBJECTDIR}/src/PPInitializeSolution_2D.o: src/PPInitializeSolution_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PPInitializeSolution_2D.o src/PPInitializeSolution_2D.c

${OBJECTDIR}/src/PPOptions_2D.o: src/PPOptions_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PPOptions_2D.o src/PPOptions_2D.c

${OBJECTDIR}/src/PPScalarFlowVariables_2D.o: src/PPScalarFlowVariables_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PPScalarFlowVariables_2D.o src/PPScalarFlowVariables_2D.c

${OBJECTDIR}/src/PPSolutionGradients_2D.o: src/PPSolutionGradients_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PPSolutionGradients_2D.o src/PPSolutionGradients_2D.c

${OBJECTDIR}/src/PPTensorFlowVariables_2D.o: src/PPTensorFlowVariables_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PPTensorFlowVariables_2D.o src/PPTensorFlowVariables_2D.c

${OBJECTDIR}/src/PPVectorFlowVariables_2D.o: src/PPVectorFlowVariables_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PPVectorFlowVariables_2D.o src/PPVectorFlowVariables_2D.c

${OBJECTDIR}/src/PostAnalysis_2D.o: src/PostAnalysis_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PostAnalysis_2D.o src/PostAnalysis_2D.c

${OBJECTDIR}/src/PostProcessingInitialize_2D.o: src/PostProcessingInitialize_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PostProcessingInitialize_2D.o src/PostProcessingInitialize_2D.c

${OBJECTDIR}/src/PostProcessing_2D.o: src/PostProcessing_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PostProcessing_2D.o src/PostProcessing_2D.c

${OBJECTDIR}/src/PrePostProcessing_2D.o: src/PrePostProcessing_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PrePostProcessing_2D.o src/PrePostProcessing_2D.c

${OBJECTDIR}/src/Read_CGNS.o: src/Read_CGNS.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Read_CGNS.o src/Read_CGNS.c

${OBJECTDIR}/src/ReferenceQuantities.o: src/ReferenceQuantities.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ReferenceQuantities.o src/ReferenceQuantities.c

${OBJECTDIR}/src/ScalarTopology_2D.o: src/ScalarTopology_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ScalarTopology_2D.o src/ScalarTopology_2D.c

${OBJECTDIR}/src/Stream_2D.o: src/Stream_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Stream_2D.o src/Stream_2D.c

${OBJECTDIR}/src/Streamline_2D.o: src/Streamline_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Streamline_2D.o src/Streamline_2D.c

${OBJECTDIR}/src/Vector.o: src/Vector.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Vector.o src/Vector.c

${OBJECTDIR}/src/VectorPlot_2D.o: src/VectorPlot_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/VectorPlot_2D.o src/VectorPlot_2D.c

${OBJECTDIR}/src/VectorTopology_2D.o: src/VectorTopology_2D.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Ilib/sources/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/VectorTopology_2D.o src/VectorTopology_2D.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/postcfd

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
