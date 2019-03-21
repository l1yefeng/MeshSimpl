#!/bin/bash
if ROOT=$(git rev-parse --show-toplevel); then
    clang-format -i -verbose -style=file \
        "${ROOT}"/src/*.{hpp,cpp} \
        "${ROOT}"/tests/test.cpp \
        "${ROOT}"/example/main.cpp

    echo
    echo "Lines longer than 90:"
    grep --color -n '.\{91\}' src/*.{hpp,cpp}
fi
