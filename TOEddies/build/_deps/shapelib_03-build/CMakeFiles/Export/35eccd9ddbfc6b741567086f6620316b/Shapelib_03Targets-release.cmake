#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Shapelib_03::shapelib_03" for configuration "Release"
set_property(TARGET Shapelib_03::shapelib_03 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Shapelib_03::shapelib_03 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libshapelib_03.a"
  )

list(APPEND _cmake_import_check_targets Shapelib_03::shapelib_03 )
list(APPEND _cmake_import_check_files_for_Shapelib_03::shapelib_03 "${_IMPORT_PREFIX}/lib/libshapelib_03.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
