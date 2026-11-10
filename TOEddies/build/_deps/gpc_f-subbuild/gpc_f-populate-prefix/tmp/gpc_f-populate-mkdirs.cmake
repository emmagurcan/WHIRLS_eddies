# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-src")
  file(MAKE_DIRECTORY "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-src")
endif()
file(MAKE_DIRECTORY
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-build"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-subbuild/gpc_f-populate-prefix"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-subbuild/gpc_f-populate-prefix/tmp"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-subbuild/gpc_f-populate-prefix/src/gpc_f-populate-stamp"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-subbuild/gpc_f-populate-prefix/src"
  "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-subbuild/gpc_f-populate-prefix/src/gpc_f-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-subbuild/gpc_f-populate-prefix/src/gpc_f-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/gpc_f-subbuild/gpc_f-populate-prefix/src/gpc_f-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
