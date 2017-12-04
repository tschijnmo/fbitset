# Download the Catch2 single header.

file(DOWNLOAD
        "https://github.com/catchorg/Catch2/releases/download/v2.0.1/catch.hpp"
        "${CMAKE_CURRENT_BINARY_DIR}/catch.hpp"
        )

# The main test driver.
add_executable(testmain
        test.cpp
        )
target_include_directories(testmain
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
        PRIVATE "${PROJECT_SOURCE_DIR}/include"
        )

