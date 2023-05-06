set(HEADER_LIST "${Numerus_SOURCE_DIR}/include/BigInt.hpp")

add_library(BigInt BigInt.cpp ${HEADER_LIST})

target_include_directories(BigInt PUBLIC ../include)

target_compile_features(BigInt PUBLIC cxx_std_20)
target_link_libraries(BigInt armadillo)

add_executable(driver driver.cpp)

target_link_libraries(driver BigInt)

set_target_properties(driver PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(CMAKE_CXX_STANDARD 14)

source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
