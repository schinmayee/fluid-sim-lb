#!/usr/bin/env bash

pkill canary
cd ${HOME}/research/runtimes/canary
rm -f controller.txt
rm -f /tmp/canary_controller*
./build/src/canary_controller --v=1 &> /dev/null
