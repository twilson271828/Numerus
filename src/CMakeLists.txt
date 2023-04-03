set(HEADER_LIST "${Numerus_SOURCE_DIR}/include/BigInt.hpp")



add_library(BigInt BigInt.cpp ${HEADER_LIST})

target_include_directories(BigInt PUBLIC ../include)

target_compile_features(BigInt PUBLIC cxx_std_11)
target_link_libraries(BigInt)

add_executable(driver driver.cpp)

target_link_libraries(driver BigInt)
set(CMAKE_CXX_STANDARD 14)

enable_testing()


#message("Supported features = ${CMAKE_CXX_COMPILE_FEATURES}")


source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
