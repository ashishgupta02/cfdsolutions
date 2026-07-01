# EOS library
file(GLOB EOS_C_SOURCES src/*.c)
file(GLOB EOS_CPP_SOURCES src/*.cpp)
add_library(EOS ${EOS_C_SOURCES} ${EOS_CPP_SOURCES})

target_include_directories(EOS PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/MATH/include
  ${CMAKE_SOURCE_DIR}/UTILS/include
  ${CMAKE_SOURCE_DIR}/NISTThermo/NISTThermo/include
)

target_link_libraries(EOS PUBLIC MATH UTILS NISTThermo)

target_compile_features(EOS PUBLIC c_std_99 cxx_std_11)

# Apply project-wide compile options
target_link_libraries(EOS PUBLIC project_compile_options)

if(ENABLE_MPI)
  target_compile_definitions(EOS PRIVATE HAVE_MPI)
endif()

set_target_properties(EOS PROPERTIES POSITION_INDEPENDENT_CODE ON)
