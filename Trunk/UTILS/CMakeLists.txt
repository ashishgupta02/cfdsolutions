# UTILS library
file(GLOB UTILS_C_SOURCES src/*.c)
file(GLOB UTILS_CPP_SOURCES src/*.cpp)
add_library(UTILS ${UTILS_C_SOURCES} ${UTILS_CPP_SOURCES})

target_include_directories(UTILS PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}
)

target_compile_features(UTILS PUBLIC c_std_99 cxx_std_11)

# Apply project-wide compile options
target_link_libraries(UTILS PUBLIC project_compile_options)

if(ENABLE_MPI)
  target_compile_definitions(UTILS PRIVATE HAVE_MPI)
endif()

set_target_properties(UTILS PROPERTIES POSITION_INDEPENDENT_CODE ON)
