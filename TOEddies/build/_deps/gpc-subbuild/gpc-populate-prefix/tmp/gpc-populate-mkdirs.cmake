# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-src")
  file(MAKE_DIRECTORY "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-src")
endif()
file(MAKE_DIRECTORY
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-build"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-subbuild/gpc-populate-prefix"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-subbuild/gpc-populate-prefix/tmp"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-subbuild/gpc-populate-prefix/src/gpc-populate-stamp"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-subbuild/gpc-populate-prefix/src"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-subbuild/gpc-populate-prefix/src/gpc-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-subbuild/gpc-populate-prefix/src/gpc-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc-subbuild/gpc-populate-prefix/src/gpc-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
