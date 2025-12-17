cmake_minimum_required(VERSION 3.20)

project(lncrna_mers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

set(seqan3_git_tag "3.0.3")
message(STATUS "Fetching SeqAn3 ${seqan3_git_tag}")

FetchContent_Declare(
    seqan3_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
    GIT_TAG        "${seqan3_git_tag}"
)

FetchContent_MakeAvailable(seqan3_fetch_content)

add_executable(lncrna_mers
    src/main.cpp
    src/logger.cpp
    src/fragmenter.cpp
    src/ibf_index.cpp
    src/query_processor.cpp
)

target_include_directories(lncrna_mers PRIVATE src)

target_link_libraries(lncrna_mers
    PUBLIC
        seqan3::seqan3
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(lncrna_mers PRIVATE -Wall -Wextra -Wpedantic -O3)
endif()
