#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Contour_531::contour_531" for configuration "Release"
set_property(TARGET Contour_531::contour_531 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Contour_531::contour_531 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcontour_531.a"
  )

list(APPEND _cmake_import_check_targets Contour_531::contour_531 )
list(APPEND _cmake_import_check_files_for_Contour_531::contour_531 "${_IMPORT_PREFIX}/lib/libcontour_531.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
