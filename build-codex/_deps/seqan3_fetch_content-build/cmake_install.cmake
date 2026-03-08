# Install script for directory: /Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/seqan3" TYPE FILE FILES
    "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/CHANGELOG.md"
    "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/CODE_OF_CONDUCT.md"
    "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/CONTRIBUTING.md"
    "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/LICENSE.md"
    "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/README.md"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/seqan3" TYPE FILE FILES
    "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/build_system/seqan3-config.cmake"
    "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/build_system/seqan3-config-version.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/include/seqan3")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/seqan3/submodules/sdsl-lite" TYPE DIRECTORY FILES "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/submodules/sdsl-lite/include")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/seqan3/submodules/range-v3" TYPE DIRECTORY FILES "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/submodules/range-v3/include")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/seqan3/submodules/cereal" TYPE DIRECTORY FILES "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build/_deps/seqan3_fetch_content-src/submodules/cereal/include")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/build-codex/_deps/seqan3_fetch_content-build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
