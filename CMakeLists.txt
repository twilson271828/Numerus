cmake_minimum_required(VERSION 3.16)

project(
    Numerus
    VERSION 1.0
    LANGUAGES CXX)

    enable_testing()

    add_subdirectory(src)

    add_subdirectory(include)

    add_subdirectory(tests)

    


    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
    include(${CMAKE_MODULE_PATH}/CodeCoverage.cmake)

    
    if(CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    endif()
    
    # Find gcovr
    find_program(GCOVR gcovr)
    if(NOT GCOVR)
        message(FATAL_ERROR "gcovr not found! Please install it.")
    endif()
    
    add_custom_target(coverage
    COMMAND ${CMAKE_COMMAND} -E env GCOV_PREFIX=${CMAKE_BINARY_DIR}/CMakeFiles/cryptome.dir
    COMMAND ${CMAKE_COMMAND} -E env GCOV_PREFIX_STRIP=3
    COMMAND ${CMAKE_COMMAND} -E env GCOVR_CONFIG=${CMAKE_SOURCE_DIR}/.gcovr.cfg
    COMMAND gcovr -r ${CMAKE_SOURCE_DIR} --html --html-details -o coverage.html
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS cryptome
)

