# Install script for directory: /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-src

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/Attributes/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/Datasets/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/Groups/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/Variables/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/modules/netcdf95.mod")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/libnetcdf95.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnetcdf95.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnetcdf95.a")
    execute_process(COMMAND "/Users/emmagurcan/miniforge3/envs/toeddies_pkgs/bin/arm64-apple-darwin20.0.0-ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnetcdf95.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/NetCDF95" TYPE FILE FILES
    "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/NetCDF95Config.cmake"
    "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/FindNetCDF_Fortran.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/NetCDF95/NetCDF95Targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/NetCDF95/NetCDF95Targets.cmake"
         "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/CMakeFiles/Export/43f992f9e44a96cbe9eefd02446eec55/NetCDF95Targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/NetCDF95/NetCDF95Targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/NetCDF95/NetCDF95Targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/NetCDF95" TYPE FILE FILES "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/CMakeFiles/Export/43f992f9e44a96cbe9eefd02446eec55/NetCDF95Targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/NetCDF95" TYPE FILE FILES "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/CMakeFiles/Export/43f992f9e44a96cbe9eefd02446eec55/NetCDF95Targets-release.cmake")
  endif()
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/_deps/netcdf95-build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
