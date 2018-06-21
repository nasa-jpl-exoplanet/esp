#! /usr/bin/env bash

. .ci/util.sh

state="pending" # "success" "pending" "failure" "error"
description="Compliancy check to ensure DAWGIE is satisfied"
context="continuous-integration/03/esp-dawgie-compliance"

post_state "$context" "$description" "$state"

if current_state
then
    python3 -m dawgie.tools.compliant
    state=`get_state`
fi

post_state "$context" "$description" "$state"
