# MESHQA executable

if(NOT ENABLE_CGNS)
  message(STATUS "Skipping MESHQA because ENABLE_CGNS is OFF")
  return()
endif()

file(GLOB MESHQA_SOURCES src/*.cpp)
add_executable(meshqa ${MESHQA_SOURCES})

target_include_directories(meshqa PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/CFDDB/include
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}
)

target_link_libraries(meshqa PRIVATE CFDDB UTILS)

target_compile_features(meshqa PRIVATE cxx_std_11)
target_compile_features(meshqa PUBLIC c_std_99 cxx_std_11)

# Apply project-wide compile options
target_link_libraries(meshqa PRIVATE project_compile_options)

target_include_directories(meshqa PRIVATE ${CGNS_INCLUDE_DIR})

target_link_libraries(meshqa PRIVATE ${CGNS_LIBRARY})
