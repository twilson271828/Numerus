cmake_minimum_required(VERSION 3.16)

project(
    Numerus
    VERSION 1.0
    LANGUAGES CXX)

    enable_testing()

    add_subdirectory(src)
    add_subdirectory(include)
    add_subdirectory(tests)
    set(CMAKE_VERBOSE_MAKEFILE ON)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
    include(${CMAKE_MODULE_PATH}/CodeCoverage.cmake)

    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

    endif()

    add_test(NAME BigIntUnitTest COMMAND BigIntUnitTest)
    set_target_properties(BigIntUnitTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testbin")

    #message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    #message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

    
    # Find gcovr
    find_program(GCOVR gcovr)
    if(NOT GCOVR)
        message(FATAL_ERROR "gcovr not found! Please install it.")
    endif()
    
#message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
