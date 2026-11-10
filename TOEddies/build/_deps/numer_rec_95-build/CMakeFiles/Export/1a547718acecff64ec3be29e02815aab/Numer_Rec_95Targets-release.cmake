#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Numer_Rec_95::numer_rec_95" for configuration "Release"
set_property(TARGET Numer_Rec_95::numer_rec_95 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Numer_Rec_95::numer_rec_95 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libnumer_rec_95.a"
  )

list(APPEND _cmake_import_check_targets Numer_Rec_95::numer_rec_95 )
list(APPEND _cmake_import_check_files_for_Numer_Rec_95::numer_rec_95 "${_IMPORT_PREFIX}/lib/libnumer_rec_95.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
