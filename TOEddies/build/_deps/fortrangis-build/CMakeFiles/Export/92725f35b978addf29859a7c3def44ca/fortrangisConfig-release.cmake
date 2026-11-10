#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fortrangis::fortranc" for configuration "Release"
set_property(TARGET fortrangis::fortranc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(fortrangis::fortranc PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libfortranc.a"
  )

list(APPEND _cmake_import_check_targets fortrangis::fortranc )
list(APPEND _cmake_import_check_files_for_fortrangis::fortranc "${_IMPORT_PREFIX}/lib/libfortranc.a" )

# Import target "fortrangis::fortranGIS" for configuration "Release"
set_property(TARGET fortrangis::fortranGIS APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(fortrangis::fortranGIS PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libfortranGIS.a"
  )

list(APPEND _cmake_import_check_targets fortrangis::fortranGIS )
list(APPEND _cmake_import_check_files_for_fortrangis::fortranGIS "${_IMPORT_PREFIX}/lib/libfortranGIS.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
