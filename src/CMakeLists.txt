
SET(SOURCES aira.cpp cluster.cpp cluster.h env.h errors.h lockingvars.h problem.cpp problem.h result.cpp result.h sense.h solutions.cpp solutions.h symgroup.cpp symgroup_extern.h symgroup.h thread.cpp thread.h types.h)

ADD_EXECUTABLE(moip_aira ${SOURCES})
TARGET_LINK_LIBRARIES(moip_aira ${Boost_PROGRAM_OPTIONS_LIBRARY} "cplex")
