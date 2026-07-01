# MPFLOW3D executable
file(GLOB MPFLOW3D_C_SOURCES src/*.c)
file(GLOB MPFLOW3D_CPP_SOURCES src/*.cpp)
add_executable(mpflow3d main.cpp ${MPFLOW3D_C_SOURCES} ${MPFLOW3D_CPP_SOURCES})

target_include_directories(mpflow3d PRIVATE
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}/MATH/include
  ${CMAKE_SOURCE_DIR}/EOS/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(mpflow3d PRIVATE
  UTILS
  MATH
  EOS
  NISTThermo
)

# Apply project-wide compile options
target_link_libraries(mpflow3d PRIVATE project_compile_options)

if(ENABLE_MPI)
  target_compile_definitions(mpflow3d PRIVATE HAVE_MPI)
  target_link_libraries(mpflow3d PRIVATE MPI::MPI_C MPI::MPI_CXX)
endif()
