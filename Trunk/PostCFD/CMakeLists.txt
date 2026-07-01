# PostCFD executable

if(NOT ENABLE_CGNS)
  message(STATUS "Skipping PostCFD because ENABLE_CGNS is OFF")
  return()
endif()

file(GLOB POSTCFD_SOURCES src/*.c)
add_executable(postcfd ${POSTCFD_SOURCES})

target_include_directories(postcfd PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/sources/include
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/UTILS/include
)

target_link_libraries(postcfd PRIVATE UTILS ${CGNS_LIBRARY})

target_compile_features(postcfd PRIVATE c_std_99)
target_compile_features(postcfd PRIVATE cxx_std_11)

# Apply project-wide compile options
target_link_libraries(postcfd PRIVATE project_compile_options)
