# Taken from the KDE4 cmake modules directory as FindGMP.cmake, then editted
# William Pettersson, 2017.

# Try to find the CPLEX librairies
#  CPLEX_FOUND - system has CPLEX
#  CPLEX_INCLUDE_DIR - the CPLEX include directory
#  CPLEX_LIBRARY - Libraries needed to use CPLEX
#  CPLEX_LIBRARY_PATH - Path to CPLEX libraries
#
# To tell this module where to search (i.e. non-standard locations) use:
#
# CPLEX_ROOT = /some/other/path/IBM/ILOG/CPLEX_Studio127

# Copyright (c) 2006, Laurent Montel, <montel@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


if (CPLEX_INCLUDE_DIR AND CPLEX_LIBRARIES)
  # Already in cache, be silent
  set(CPLEX_FIND_QUIETLY TRUE)
endif (CPLEX_INCLUDE_DIR AND CPLEX_LIBRARIES)

find_path(CPLEX_INCLUDE_DIR NAMES ilcplex/cplex.h
  HINTS "${CPLEX_ROOT}/cplex/include")


IF(CPLEX_INCLUDE_DIR AND EXISTS "${CPLEX_INCLUDE_DIR}/ilcplex/cplex.h")
  # TODO Detect more systems/architectures
  IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
      SET(_CPLEX_LIBRARY_HINT "x86-64_linux")
    ELSE()
      SET(_CPLEX_LIBRARY_HINT "x86_linux")
    ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

  find_library(_CPLEX_LIBRARIES NAMES cplex
    HINTS "${CPLEX_ROOT}/cplex/lib/${_CPLEX_LIBRARY_HINT}/static_pic/")

  IF(_CPLEX_LIBRARIES)
    GET_FILENAME_COMPONENT(CPLEX_LIBRARY_PATH ${_CPLEX_LIBRARIES} PATH)
    IF(CPLEX_LIBRARY_PATH)
      SET(CPLEX_LIBRARY "cplex")
      SET(CPLEX_Found)
    ENDIF(CPLEX_LIBRARY_PATH)
  ENDIF(_CPLEX_LIBRARIES)
ENDIF(CPLEX_INCLUDE_DIR AND EXISTS "${CPLEX_INCLUDE_DIR}/ilcplex/cplex.h")


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG CPLEX_INCLUDE_DIR CPLEX_LIBRARY)

mark_as_advanced(CPLEX_INCLUDE_DIR CPLEX_LIBRARY)
