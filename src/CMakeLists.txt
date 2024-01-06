set(HEADER_LIST "${Numerus_SOURCE_DIR}/include/BigInt.hpp ${Numerus_SOURCE_DIR}/include/Numerus.hpp")


add_library(BigInt BigInt.cpp ${HEADER_LIST})
add_library(Numerus Numerus.cpp ${HEADER_LIST})

target_include_directories(BigInt Numerus PUBLIC ../include)


target_compile_features(BigInt Numerus PUBLIC cxx_std_20)
target_link_libraries(BigInt armadillo)

add_executable(driver driver.cpp)
add_executable(cryptome cryptome.cpp)

#target_link_libraries(driver BigInt)

target_link_libraries(cryptome BigInt Numerus)

set_target_properties(driver cryptome PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")



source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
