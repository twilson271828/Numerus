enable_testing()

find_package(GTest REQUIRED)

add_executable(BigIntUnitTest BigIntUnitTest.cpp)

target_link_libraries(BigIntUnitTest
 PRIVATE
  GTest::GTest
  BigInt)

  add_test(BigIntUnitTest BigIntUnitTest)
  
  set_target_properties(BigIntUnitTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testbin")



