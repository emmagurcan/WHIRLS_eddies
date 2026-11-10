#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Jumble::jumble" for configuration "Release"
set_property(TARGET Jumble::jumble APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Jumble::jumble PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libjumble.a"
  )

list(APPEND _cmake_import_check_targets Jumble::jumble )
list(APPEND _cmake_import_check_files_for_Jumble::jumble "${_IMPORT_PREFIX}/lib/libjumble.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
