# Euler2D executable

file(GLOB_RECURSE EULER2D_SOLVER_SOURCES SOLVER/src/*.cpp)
file(GLOB_RECURSE EULER2D_MESH_SOURCES MESH/src/*.cpp)
file(GLOB_RECURSE EULER2D_DESIGN_SOURCES DESIGN/src/*.cpp)
# MESH/src/Mesh.cpp is a superseded 2010 standalone implementation replaced by
# Euler2D_Mesh.cpp; it is not part of the executable.
list(FILTER EULER2D_MESH_SOURCES EXCLUDE REGEX "MESH/src/Mesh\\.cpp$")
add_executable(euler2d main.cpp ${EULER2D_SOLVER_SOURCES} ${EULER2D_MESH_SOURCES} ${EULER2D_DESIGN_SOURCES})

target_include_directories(euler2d PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/SOLVER/include
  ${CMAKE_CURRENT_SOURCE_DIR}/MESH/include
  ${CMAKE_CURRENT_SOURCE_DIR}/DESIGN/include
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}/MATH/include
  ${CMAKE_SOURCE_DIR}
)

target_link_libraries(euler2d PRIVATE UTILS MATH)

# Apply project-wide compile options
target_link_libraries(euler2d PRIVATE project_compile_options)

target_compile_features(euler2d PRIVATE cxx_std_11)
