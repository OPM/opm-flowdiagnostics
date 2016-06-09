#!/bin/bash

source `dirname $0`/build-opm-flowdiagnostics.sh

OPM_COMMON_REVISION=master

build_opm_flowdiagnostics
test $? -eq 0 || exit 1

cp serial/build-opm-flowdiagnostics/testoutput.xml .
