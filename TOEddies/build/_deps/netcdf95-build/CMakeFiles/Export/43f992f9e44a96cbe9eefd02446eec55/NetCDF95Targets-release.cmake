#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "NetCDF95::netcdf95" for configuration "Release"
set_property(TARGET NetCDF95::netcdf95 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(NetCDF95::netcdf95 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libnetcdf95.a"
  )

list(APPEND _cmake_import_check_targets NetCDF95::netcdf95 )
list(APPEND _cmake_import_check_files_for_NetCDF95::netcdf95 "${_IMPORT_PREFIX}/lib/libnetcdf95.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
