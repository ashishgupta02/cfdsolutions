# UGRID2SU2 executable

add_executable(ugrid2su2 main.cpp)

target_include_directories(ugrid2su2 PRIVATE
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}
)

target_link_libraries(ugrid2su2 PRIVATE UTILS)

target_compile_features(ugrid2su2 PRIVATE cxx_std_11)
 
# Apply project-wide compile options
target_link_libraries(ugrid2su2 PRIVATE project_compile_options)
