set(HEADER_LIST "${Numerus_SOURCE_DIR}/include/BigInt.hpp")

add_library(BigInt BigInt.cpp ${HEADER_LIST})

target_include_directories(BigInt PUBLIC ../include)

target_compile_features(BigInt PUBLIC cxx_std_11)
target_link_libraries(BigInt)

add_executable(driver driver.cpp)

set_target_properties(driver PROPERTIES FOLDER ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


target_link_libraries(driver BigInt)
set(CMAKE_CXX_STANDARD 14)

source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
