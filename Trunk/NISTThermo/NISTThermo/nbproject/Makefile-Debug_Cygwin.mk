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
CND_PLATFORM=Cygwin_4.x-Windows
CND_DLIB_EXT=dll
CND_CONF=Debug_Cygwin
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/cmns.o \
	${OBJECTDIR}/src/core_anc.o \
	${OBJECTDIR}/src/core_bwr.o \
	${OBJECTDIR}/src/core_cpp.o \
	${OBJECTDIR}/src/core_de.o \
	${OBJECTDIR}/src/core_ecs.o \
	${OBJECTDIR}/src/core_feq.o \
	${OBJECTDIR}/src/core_mlt.o \
	${OBJECTDIR}/src/core_ph0.o \
	${OBJECTDIR}/src/core_pr.o \
	${OBJECTDIR}/src/core_qui.o \
	${OBJECTDIR}/src/core_stn.o \
	${OBJECTDIR}/src/flash2.o \
	${OBJECTDIR}/src/flsh_sub.o \
	${OBJECTDIR}/src/idealgas.o \
	${OBJECTDIR}/src/mix_aga8.o \
	${OBJECTDIR}/src/mix_hmx.o \
	${OBJECTDIR}/src/nist_extension_mp_pd.o \
	${OBJECTDIR}/src/nist_extension_mp_td.o \
	${OBJECTDIR}/src/nist_extension_pd.o \
	${OBJECTDIR}/src/nist_extension_td.o \
	${OBJECTDIR}/src/prop_sub.o \
	${OBJECTDIR}/src/realgas.o \
	${OBJECTDIR}/src/sat_sub.o \
	${OBJECTDIR}/src/setup.o \
	${OBJECTDIR}/src/setup2.o \
	${OBJECTDIR}/src/trns_ecs.o \
	${OBJECTDIR}/src/trns_tcx.o \
	${OBJECTDIR}/src/trns_vis.o \
	${OBJECTDIR}/src/trnsp.o \
	${OBJECTDIR}/src/utility.o


# C Compiler Flags
CFLAGS=-march=native -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2

# CC Compiler Flags
CCFLAGS=-march=native -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2
CXXFLAGS=-march=native -Wno-write-strings -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -fpermissive -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2

# Fortran Compiler Flags
FFLAGS=-march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libNISTThermo.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libNISTThermo.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libNISTThermo.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -march=native -shared

${OBJECTDIR}/src/cmns.o: src/cmns.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/cmns.o src/cmns.f

${OBJECTDIR}/src/core_anc.o: src/core_anc.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_anc.o src/core_anc.f

${OBJECTDIR}/src/core_bwr.o: src/core_bwr.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_bwr.o src/core_bwr.f

${OBJECTDIR}/src/core_cpp.o: src/core_cpp.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_cpp.o src/core_cpp.f

${OBJECTDIR}/src/core_de.o: src/core_de.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_de.o src/core_de.f

${OBJECTDIR}/src/core_ecs.o: src/core_ecs.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_ecs.o src/core_ecs.f

${OBJECTDIR}/src/core_feq.o: src/core_feq.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_feq.o src/core_feq.f

${OBJECTDIR}/src/core_mlt.o: src/core_mlt.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_mlt.o src/core_mlt.f

${OBJECTDIR}/src/core_ph0.o: src/core_ph0.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_ph0.o src/core_ph0.f

${OBJECTDIR}/src/core_pr.o: src/core_pr.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_pr.o src/core_pr.f

${OBJECTDIR}/src/core_qui.o: src/core_qui.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_qui.o src/core_qui.f

${OBJECTDIR}/src/core_stn.o: src/core_stn.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/core_stn.o src/core_stn.f

${OBJECTDIR}/src/flash2.o: src/flash2.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/flash2.o src/flash2.f

${OBJECTDIR}/src/flsh_sub.o: src/flsh_sub.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/flsh_sub.o src/flsh_sub.f

${OBJECTDIR}/src/idealgas.o: src/idealgas.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/idealgas.o src/idealgas.f

${OBJECTDIR}/src/mix_aga8.o: src/mix_aga8.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/mix_aga8.o src/mix_aga8.f

${OBJECTDIR}/src/mix_hmx.o: src/mix_hmx.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/mix_hmx.o src/mix_hmx.f

${OBJECTDIR}/src/nist_extension_mp_pd.o: src/nist_extension_mp_pd.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/nist_extension_mp_pd.o src/nist_extension_mp_pd.f

${OBJECTDIR}/src/nist_extension_mp_td.o: src/nist_extension_mp_td.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/nist_extension_mp_td.o src/nist_extension_mp_td.f

${OBJECTDIR}/src/nist_extension_pd.o: src/nist_extension_pd.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/nist_extension_pd.o src/nist_extension_pd.f

${OBJECTDIR}/src/nist_extension_td.o: src/nist_extension_td.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/nist_extension_td.o src/nist_extension_td.f

${OBJECTDIR}/src/prop_sub.o: src/prop_sub.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/prop_sub.o src/prop_sub.f

${OBJECTDIR}/src/realgas.o: src/realgas.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/realgas.o src/realgas.f

${OBJECTDIR}/src/sat_sub.o: src/sat_sub.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/sat_sub.o src/sat_sub.f

${OBJECTDIR}/src/setup.o: src/setup.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/setup.o src/setup.f

${OBJECTDIR}/src/setup2.o: src/setup2.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/setup2.o src/setup2.f

${OBJECTDIR}/src/trns_ecs.o: src/trns_ecs.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/trns_ecs.o src/trns_ecs.f

${OBJECTDIR}/src/trns_tcx.o: src/trns_tcx.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/trns_tcx.o src/trns_tcx.f

${OBJECTDIR}/src/trns_vis.o: src/trns_vis.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/trns_vis.o src/trns_vis.f

${OBJECTDIR}/src/trnsp.o: src/trnsp.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/trnsp.o src/trnsp.f

${OBJECTDIR}/src/utility.o: src/utility.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -g -Wall -march=native -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000 -Winline -mfpmath=sse -msse2 -fdefault-double-8 -fdefault-real-8 -I./include  -o ${OBJECTDIR}/src/utility.o src/utility.f

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libNISTThermo.${CND_DLIB_EXT}
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
