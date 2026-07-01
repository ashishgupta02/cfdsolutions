# Centralized compiler and configuration flags

# Set a default build type for single-config generators
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

add_library(project_compile_options INTERFACE)

# Debug/Release compile definitions
target_compile_definitions(project_compile_options INTERFACE
  $<$<CONFIG:Debug>:DEBUG>
)

# Optionally tune generated code for the building machine's CPU (matches the
# original build). Turn OFF for portable/redistributable binaries.
option(ENABLE_NATIVE_ARCH "Optimize for the building machine's CPU (-march=native)" ON)
set(_native_arch)
if(ENABLE_NATIVE_ARCH)
  set(_native_arch -march=native)
endif()

# GNU (gcc/g++) Release optimization and floating-point profile carried over
# from the original NetBeans CFLAGS/CXXFLAGS.
# WARNING: -ffast-math implies -ffinite-math-only, which lets the compiler
# assume values are never NaN/Inf. isnan()-based blow-up checks may therefore be
# optimized away in Release builds; keep this in mind when diagnosing solver
# divergence.
set(_gnu_opt
  ${_native_arch}
  -mfpmath=sse -msse2
  -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans
  -fstrict-aliasing -fomit-frame-pointer -fexpensive-optimizations
  -funroll-all-loops -ffast-math -finline-limit=150000
)

# GNU Fortran Release profile (matches FFLAGS). The NIST thermodynamics sources
# require 8-byte reals/doubles.
set(_gnu_fortran_opt
  ${_native_arch}
  -mfpmath=sse -msse2
  -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-limit=150000
)
set(_gnu_fortran_real8 -fdefault-double-8 -fdefault-real-8)

# Clang: the portable subset of the same intent (Clang rejects several of the
# GNU-only tuning switches).
set(_clang_opt ${_native_arch} -ffast-math -funroll-loops -fstrict-aliasing -fomit-frame-pointer)

target_compile_options(project_compile_options INTERFACE
  # ---- GNU ----
  "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,GNU>>:-g;-Og>"
  "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:CXX,GNU>>:-g;-Og>"
  "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-g;-Og;${_gnu_fortran_real8}>"
  "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O3;-fPIC;${_gnu_opt}>"
  "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,GNU>>:-O3;-fPIC;${_gnu_opt}>"
  "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-O3;-fPIC;${_gnu_fortran_opt};${_gnu_fortran_real8}>"

  # ---- Clang ----
  "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,Clang>>:-g;-O0>"
  "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:CXX,Clang>>:-g;-O0>"
  "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,Clang>>:-O3;-fPIC;${_clang_opt}>"
  "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,Clang>>:-O3;-fPIC;${_clang_opt}>"

  # ---- Intel (classic & oneAPI/LLVM) ----
  "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,Intel,IntelLLVM>>:-g;-O0>"
  "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:CXX,Intel,IntelLLVM>>:-g;-O0>"
  "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,Intel,IntelLLVM>>:-xHost;-O3;-fPIC>"
  "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,Intel,IntelLLVM>>:-xHost;-O3;-fPIC>"
)

# Intel Fortran: 8-byte reals/doubles
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
  target_compile_options(project_compile_options INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-r8;-real-size;64>)
endif()

# Optionally enable stricter warnings (can be toggled by users)
option(ENABLE_STRICT_WARNINGS "Enable strict compiler warnings" ON)
if(ENABLE_STRICT_WARNINGS)
  target_compile_options(project_compile_options INTERFACE
    $<$<OR:$<C_COMPILER_ID:GNU>,$<C_COMPILER_ID:Clang>>:-Wall;-Wextra;-Wpedantic>
    $<$<C_COMPILER_ID:Intel>:-warn;all>
  )
endif()
