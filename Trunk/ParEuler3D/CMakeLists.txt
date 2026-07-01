# ParEuler3D executable

if(NOT ENABLE_MPI)
  message(STATUS "Skipping ParEuler3D because ENABLE_MPI is OFF")
  return()
endif()

file(GLOB PARCELLIO_SOURCES MESHIO/src/*.cpp)
file(GLOB PARCOMM_MPI_SOURCES CommMPI/src/*.cpp)
add_executable(par_euler3d main.cpp Commons.cpp ${PARCELLIO_SOURCES} ${PARCOMM_MPI_SOURCES})

target_include_directories(par_euler3d PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/MESHIO/include
  ${CMAKE_CURRENT_SOURCE_DIR}/CommMPI/include
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}/MATH/include
  ${CMAKE_SOURCE_DIR}
)

target_link_libraries(par_euler3d PRIVATE UTILS MATH)

target_compile_features(par_euler3d PRIVATE c_std_99 cxx_std_11)

# Apply project-wide compile options
target_link_libraries(par_euler3d PRIVATE project_compile_options)

if(ENABLE_MPI)
  target_compile_definitions(par_euler3d PRIVATE HAVE_MPI)
  target_link_libraries(par_euler3d PRIVATE MPI::MPI_CXX)
endif()
