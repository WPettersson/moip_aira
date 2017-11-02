#!/usr/bin/env bash

EXECUTABLE=$1
TEST=$2
OPTS=$3
TESTNAME=$(basename ${TEST} .lp)
TESTDIR=$(dirname ${TEST})
OUTFILE=$(mktemp ${TESTNAME}.XXX)
${EXECUTABLE} -p ${TEST} -o ${OUTFILE} ${OPTS}
diff -w -I 'seconds\|solved\|Using' ${TESTDIR}/${TESTNAME}.out ${OUTFILE}
RES=$?
rm ${OUTFILE}
exit ${RES}
