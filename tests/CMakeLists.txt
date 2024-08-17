find_package(GTest REQUIRED)

add_executable(BigIntUnitTest BigIntUnitTest.cpp)

target_link_libraries(BigIntUnitTest
 PRIVATE
  GTest::GTest
  BigInt)

  
  add_test(NAME BigIntUnitTest COMMAND BigIntUnitTest)
  
target_compile_options(BigIntUnitTest PRIVATE --coverage)
target_link_libraries(BigIntUnitTest PRIVATE --coverage)

  set_target_properties(BigIntUnitTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testbin")

  

