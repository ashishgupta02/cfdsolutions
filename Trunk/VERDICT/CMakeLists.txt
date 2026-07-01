# VERDICT math utility library

file(GLOB VERDICT_SOURCES src/*.cpp)
add_library(VERDICT STATIC ${VERDICT_SOURCES})

target_include_directories(VERDICT PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}
)

target_compile_features(VERDICT PUBLIC cxx_std_11)

# Apply project-wide compile options
target_link_libraries(VERDICT PUBLIC project_compile_options)
