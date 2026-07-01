# NISTThermo library
file(GLOB NISTTHERMO_FORTRAN src/*.f)
# pass_ftn.f only defines the optional DLL entry points (RPVersion, SETUPdll,
# ...) using Lahey/Fujitsu-specific `dll_export` statements that GNU/Intel
# Fortran cannot compile. The library callers use the underlying routines
# directly, so this wrapper is not needed.
list(FILTER NISTTHERMO_FORTRAN EXCLUDE REGEX "pass_ftn\\.f$")
if(NISTTHERMO_FORTRAN AND ENABLE_FORTRAN)
  add_library(NISTThermo ${NISTTHERMO_FORTRAN})
  set_target_properties(NISTThermo PROPERTIES LINKER_LANGUAGE Fortran)

  target_include_directories(NISTThermo PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
  )

  # Apply project-wide compile options. GNU/Intel 8-byte real/double flags are
  # handled centrally in cmake/CompilerFlags.cmake.
  target_link_libraries(NISTThermo PUBLIC project_compile_options)

  # Best-effort 8-byte real/double mapping for compilers not covered centrally.
  if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NAG|XL")
    target_compile_options(NISTThermo PRIVATE -fdefault-double-8 -fdefault-real-8)
  endif()

  set_target_properties(NISTThermo PROPERTIES POSITION_INDEPENDENT_CODE ON)
else()
  message(STATUS "Skipping NISTThermo Fortran library because no Fortran compiler is available")
endif()
