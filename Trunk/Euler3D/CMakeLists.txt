# Euler3D executable

file(GLOB EULER3D_SOURCES src/*.cpp)
add_executable(euler3d main.cpp ${EULER3D_SOURCES})

target_include_directories(euler3d PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}/MATH/include
  ${CMAKE_SOURCE_DIR}
)

target_link_libraries(euler3d PRIVATE UTILS MATH)

# Apply project-wide compile options
target_link_libraries(euler3d PRIVATE project_compile_options)

target_compile_features(euler3d PRIVATE cxx_std_11)

if(ENABLE_MPI)
  target_compile_definitions(euler3d PRIVATE HAVE_MPI)
  target_link_libraries(euler3d PRIVATE MPI::MPI_CXX)
endif()
