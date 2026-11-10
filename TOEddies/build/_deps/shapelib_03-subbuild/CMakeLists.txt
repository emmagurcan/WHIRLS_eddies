# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 4.4.5)

# Reject any attempt to use a toolchain file. We must not use one because
# we could be downloading it here. If the CMAKE_TOOLCHAIN_FILE environment
# variable is set, the cache variable will have been initialized from it.
unset(CMAKE_TOOLCHAIN_FILE CACHE)
unset(ENV{CMAKE_TOOLCHAIN_FILE})

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(shapelib_03-populate NONE)


# Pass through things we've already detected in the main project to avoid
# paying the cost of redetecting them again in ExternalProject_Add()
set(GIT_EXECUTABLE [==[/Users/emmagurcan/miniforge3/envs/toeddies_pkgs/bin/git]==])
set(Git_VERSION [==[2.56.0]==])
set_property(GLOBAL PROPERTY _CMAKE_FindGit_GIT_EXECUTABLE_VERSION
  [==[/Users/emmagurcan/miniforge3/envs/toeddies_pkgs/bin/git;2.56.0]==]
)


include(ExternalProject)
ExternalProject_Add(shapelib_03-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/lguez/Shapelib_03.git" "EXTERNALPROJECT_INTERNAL_ARGUMENT_SEPARATOR" "GIT_TAG" "ca08580"
                    SOURCE_DIR          "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/shapelib_03-src"
                    BINARY_DIR          "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/shapelib_03-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


