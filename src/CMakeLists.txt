set(HEADER_LIST "${Numerus_SOURCE_DIR}/include/BigInt.hpp")

add_library(BigInt BigInt.cpp ${HEADER_LIST})

target_include_directories(BigInt PUBLIC ../include)

target_compile_features(BigInt PUBLIC cxx_std_20)
target_compile_options(BigInt PUBLIC --coverage)

target_link_libraries(BigInt gcov)


add_executable(driver driver.cpp)
add_executable(cryptome cryptome.cpp)
target_link_libraries(cryptome BigInt)
target_link_libraries(driver BigInt)

set_target_properties(driver cryptome PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
