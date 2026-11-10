#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "GPC::GPC" for configuration "Release"
set_property(TARGET GPC::GPC APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GPC::GPC PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libGPC.a"
  )

list(APPEND _cmake_import_check_targets GPC::GPC )
list(APPEND _cmake_import_check_files_for_GPC::GPC "${_IMPORT_PREFIX}/lib/libGPC.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
