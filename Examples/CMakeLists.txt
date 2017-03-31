FILE(GLOB TESTS *.lp)
FOREACH(TESTFILE ${TESTS})
  GET_FILENAME_COMPONENT(TESTNAME ${TESTFILE} NAME_WE)
  ADD_TEST(NAME ${TESTNAME} COMMAND
    "${PROJECT_SOURCE_DIR}/scripts/checkResults.sh"
    $<TARGET_FILE:aira>
    "${TESTFILE}")
  GET_FILENAME_COMPONENT(TESTNAME ${TESTFILE} NAME_WE)
  ADD_TEST(NAME "${TESTNAME}-group2" COMMAND
    "${PROJECT_SOURCE_DIR}/scripts/checkResults.sh"
    $<TARGET_FILE:aira>
    "${TESTFILE}"
    "-t 2")
  ADD_TEST(NAME "${TESTNAME}-spread2" COMMAND
    "${PROJECT_SOURCE_DIR}/scripts/checkResults.sh"
    $<TARGET_FILE:aira>
    "${TESTFILE}"
    "-t 2 -s")
  ADD_TEST(NAME "${TESTNAME}-flat2" COMMAND
    "${PROJECT_SOURCE_DIR}/scripts/checkResults.sh"
    $<TARGET_FILE:aira>
    "${TESTFILE}"
    "-t 2 --split")
  ADD_TEST(NAME "${TESTNAME}-normal2" COMMAND
    "${PROJECT_SOURCE_DIR}/scripts/checkResults.sh"
    $<TARGET_FILE:aira>
    "${TESTFILE}"
    "-t 2 --split --split-normal")
ENDFOREACH(TESTFILE)
