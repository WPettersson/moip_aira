ADD_CUSTOM_COMMAND(
  OUTPUT symgroup.cpp
  COMMAND "python" ARGS "${CMAKE_CURRENT_SOURCE_DIR}/mk_symgroup.py" "${NUM_OF_SYMGROUPS}"
)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
SET(SOURCES
  aira.cpp
  cluster.cpp
  problem.cpp
  result.cpp
  solutions.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/symgroup.cpp
  thread.cpp
  hash.cpp)



FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/mkhash.sh"
  "GIT_HASH=$(cd ${CMAKE_CURRENT_SOURCE_DIR} && git rev-parse HEAD)\necho \"#include <string>\nstd::string HASH = \\\"\${GIT_HASH}\\\";\" > \"\${1}/hash.cpp\""
  )

ADD_CUSTOM_TARGET(
  ALWAYS_RUN ALL
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/always_missing.h"
  )

ADD_CUSTOM_COMMAND(
  OUTPUT
  "${CMAKE_CURRENT_BINARY_DIR}/always_missing.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/hash.cpp"
  COMMAND "chmod" ARGS "u+x" "${CMAKE_CURRENT_BINARY_DIR}/mkhash.sh"
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/mkhash.sh" ARGS "${CMAKE_CURRENT_SOURCE_DIR}"
  )


ADD_EXECUTABLE(aira ${SOURCES})
TARGET_LINK_LIBRARIES(aira ${Boost_PROGRAM_OPTIONS_LIBRARY} ${CPLEX_LIBRARY})
