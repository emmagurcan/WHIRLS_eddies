#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Geometry::geometry" for configuration "Release"
set_property(TARGET Geometry::geometry APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Geometry::geometry PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgeometry.a"
  )

list(APPEND _cmake_import_check_targets Geometry::geometry )
list(APPEND _cmake_import_check_files_for_Geometry::geometry "${_IMPORT_PREFIX}/lib/libgeometry.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
