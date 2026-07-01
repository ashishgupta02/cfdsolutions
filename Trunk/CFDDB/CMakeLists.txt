# CFDDB library

file(GLOB CFDDB_C_SOURCES src/*.c)
file(GLOB CFDDB_CPP_SOURCES src/*.cpp)
list(FILTER CFDDB_C_SOURCES EXCLUDE REGEX "netcdfIO\\.c$")
set(CFDDB_SOURCES ${CFDDB_C_SOURCES} ${CFDDB_CPP_SOURCES})

add_library(CFDDB STATIC ${CFDDB_SOURCES})

target_include_directories(CFDDB PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}/PostCFD/lib/sources/include
)

target_link_libraries(CFDDB PUBLIC UTILS)

target_compile_features(CFDDB PUBLIC c_std_99 cxx_std_11)

# Apply project-wide compile options
target_link_libraries(CFDDB PUBLIC project_compile_options)

if(ENABLE_NETCDF)
  list(APPEND CFDDB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdfIO.c)
  target_sources(CFDDB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdfIO.c)
  target_include_directories(CFDDB PRIVATE ${NETCDF_INCLUDE_DIR})
  target_link_libraries(CFDDB PRIVATE ${NETCDF_LIBRARY})
  target_compile_definitions(CFDDB PUBLIC ENABLE_NETCDF)
endif()

if(ENABLE_CGNS)
  target_include_directories(CFDDB PRIVATE ${CGNS_INCLUDE_DIR})
  target_link_libraries(CFDDB PRIVATE ${CGNS_LIBRARY})
  target_compile_definitions(CFDDB PUBLIC ENABLE_CGNS)
endif()
