cmake_minimum_required(VERSION 3.16)

project(
    Numerus
    VERSION 1.0
    LANGUAGES CXX)

    add_subdirectory(src)

    add_subdirectory(include)

    add_subdirectory(tests)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

    
