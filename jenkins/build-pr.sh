#!/bin/bash

source `dirname $0`/build-opm-flowdiagnostics.sh

# Upstream revisions
OPM_COMMON_REVISION=master

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  OPM_COMMON_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

echo "Building with opm-common=$OPM_COMMON_REVISION opm-flowdiagnostics=$sha1"

build_opm_flowdiagnostics
test $? -eq 0 || exit 1

# If no downstream builds we are done
if ! grep -q "with downstreams" <<< $ghprbCommentBody
then
  cp serial/build-opm-flowdiagnostics/testoutput.xml .
  exit 0
fi

# Downstream revisions
declare -a downstreams
downstreams=(opm-parser
             opm-material
             opm-core
             opm-flowdiagnostics-applications)

declare -A downstreamRev
downstreamRev[opm-parser]=master
downstreamRev[opm-material]=master
downstreamRev[opm-core]=master
downstreamRev[opm-flowdiagnostics-applications]=master

# Build ERT
ERT_REVISION=master
if grep -q "ert=" <<< $ghprbCommentBody
then
  ERT_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ert=([0-9]+).*/\1/g'`/merge
fi

pushd .
mkdir -p $WORKSPACE/deps/ert
cd $WORKSPACE/deps/ert
git init .
git remote add origin https://github.com/Ensembles/ert
git fetch --depth 1 origin $ERT_REVISION:branch_to_build
test $? -eq 0 || exit 1
git checkout branch_to_build
popd

pushd .
mkdir -p serial/build-ert
cd serial/build-ert
cmake $WORKSPACE/deps/ert/devel -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install
cmake --build . --target install
test $? -eq 0 || exit 1
popd

build_downstreams opm-flowdiagnostics

test $? -eq 0 || exit 1
