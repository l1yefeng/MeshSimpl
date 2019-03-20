#!/bin/bash
if ROOT=$(git rev-parse --show-toplevel); then
    clang-format -i -verbose -style=file \
        "${ROOT}"/src/*.{hpp,cpp} \
        "${ROOT}"/tests/test.cpp \
        "${ROOT}"/example/main.cpp
fi
