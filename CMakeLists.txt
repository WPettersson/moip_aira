CMAKE_MINIMUM_REQUIRED (VERSION 3.4)
PROJECT (moip_aira)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")

# We need to determine threading flags *before* testing for C++11 std::thread.
FIND_PACKAGE(Threads REQUIRED)
IF(NOT Threads_FOUND)
  MESSAGE(FATAL_ERROR "moip_aira requires threading support.")
ENDIF(NOT Threads_FOUND)
IF(CMAKE_COMPILER_IS_GNUCXX)
  # The FindThreads module picks up -lpthread, not -pthread.
  # This is unfortunate, because std::thread must be built with -pthread
  # under g++; otherwise it builds but throws an exception at runtime.
  SET(CMAKE_CXX_FLAGS "-pthread ${CMAKE_CXX_FLAGS}")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

INCLUDE(CheckCXX11Features)
SET(HAVE_CXX11_FEATURES TRUE)
NEED_CXX11_FEATURE("auto")
NEED_CXX11_FEATURE("constexpr")
NEED_CXX11_FEATURE("nullptr")

IF(NOT HAVE_CXX11_FEATURES)
  MESSAGE(FATAL_ERROR "moip_aira requires some C++11 support to compile.")
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_SWITCH}")

FIND_PACKAGE(CPLEX REQUIRED)

INCLUDE_DIRECTORIES(${CPLEX_INCLUDE_DIR})
LINK_DIRECTORIES(${CPLEX_LIBRARY_PATH})

FIND_PACKAGE(Boost 1.57.0 COMPONENTS program_options)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

ADD_SUBDIRECTORY(src)

IF(TESTSUITE)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/Examples)
ENDIF(TESTSUITE)
