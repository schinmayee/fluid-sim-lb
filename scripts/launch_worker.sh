#!/usr/bin/env bash

cd ${HOME}/research/runtimes/canary
rm -f worker_log.txt
rm -f /tmp/canary_worker*
ulimit -c unlimited
./build/src/canary_worker --controller_host $1 --worker_threads $2 &>> worker_log.txt
