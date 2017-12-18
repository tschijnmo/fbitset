# Download the Catch2 single header.

file(DOWNLOAD
        "https://github.com/catchorg/Catch2/releases/download/v2.0.1/catch.hpp"
        "${CMAKE_CURRENT_BINARY_DIR}/catch.hpp"
        )


# The main test driver.
add_executable(testmain
        testmain.cpp
        basic.cpp
        )

target_include_directories(testmain
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
        )

# Tests supposed to crash with assertion failure.
add_executable(crash_no_ext crash_no_ext.cpp)

# Tests added for the sake of peeking at the assembly output from the compiler.
add_library(testopt testopt.cpp)


# The benchmark for looping over indices.
add_library(testloop testloop.cpp)
target_include_directories(testloop
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
        )
add_executable(testloopbits testloopbits.cpp)
target_link_libraries(testloopbits testloop)
add_executable(testloopvector testloopvector.cpp)
target_link_libraries(testloopvector testloop)

