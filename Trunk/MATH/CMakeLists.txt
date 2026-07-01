# MATH library
file(GLOB MATH_C_SOURCES src/*.c)
file(GLOB MATH_CPP_SOURCES src/*.cpp)
add_library(MATH ${MATH_C_SOURCES} ${MATH_CPP_SOURCES})

target_include_directories(MATH PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/UTILS/include
)

target_link_libraries(MATH PUBLIC UTILS)

target_compile_features(MATH PUBLIC c_std_99 cxx_std_11)

# Apply project-wide compile options
target_link_libraries(MATH PUBLIC project_compile_options)

if(ENABLE_MPI)
  target_compile_definitions(MATH PRIVATE HAVE_MPI)
endif()

set_target_properties(MATH PROPERTIES POSITION_INDEPENDENT_CODE ON)
