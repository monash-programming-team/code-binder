#!/bin/bash

compile() {
    g++ ${1}.cpp -std=c++14 -coverage -Wall -Wextra -Werror -o ${1} || exit 1
}

EXIT_STATUS=0

for X in *-tests.cpp; do
    FILE=$(basename $X .cpp)
    SHORT=$(basename $FILE -tests)
    compile ${FILE}
    ./${FILE} -x 10 -d yes || EXIT_STATUS=1
    gcov -n -o . ${FILE} > /dev/null
done

exit ${EXIT_STATUS}
