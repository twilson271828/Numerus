set(HEADER_LIST "${Numerus_SOURCE_DIR}/include/BigInt.hpp" "${Numerus_SOURCE_DIR}/include/maths.h")

add_library(BigInt SHARED BigInt.cpp ${HEADER_LIST})

target_include_directories(BigInt PUBLIC ../include)

#Set the output directory for the BigInt library
set_target_properties(BigInt PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


target_compile_features(BigInt PUBLIC cxx_std_20)
target_compile_options(BigInt PUBLIC --coverage)

target_link_libraries(BigInt gcov)


add_executable(cryptome cryptome.cpp)
add_executable(prime_boss prime_boss.cpp)

target_link_libraries(cryptome BigInt)


# Link libraries: ZLIB + OpenMP
target_link_libraries(prime_boss
    PRIVATE
    ZLIB::ZLIB
    OpenMP::OpenMP_CXX
    gcov
    BigInt
)


# (Optional) More strict compilation flags

target_compile_options(prime_boss PRIVATE -Wall -Wextra -O3)

set_target_properties(cryptome prime_boss PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
