find_package(GTest REQUIRED)
if (NOT GTEST_LIBRARY)
  message("GTest not found!")
endif()

add_executable(BigIntUnitTest BigIntUnitTest.cpp)

target_link_libraries(BigIntUnitTest
 PRIVATE
  GTest::GTest
  BigInt)

  add_test(BigIntUnitTest BigIntUnitTest)
  
  set_target_properties(BigIntUnitTest PROPERTIES FOLDER ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testbin)



