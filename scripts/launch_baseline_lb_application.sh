#!/usr/bin/env bash

cd ${HOME}/research/runtimes/canary

echo "Launchiing application"

params="
  --launch_application=${HOME}/research/fluid-simulations/projects/flip-current-lb/build/water_app.so \
  --launch_placement_algorithm=${PLACEMENT} --launch_num_worker=${NUM_WORKERS} \
  total_workers=${NUM_WORKERS} \
  init=${INIT} frames=${FRAMES} ngrid=${NGRID} \
  partition_x=${X} partition_y=${Y} partition_z=${Z} \
  ax=${AX} ay=${AY} az=${AZ} \
  current_load_factor=${F} \
  processors=${P} horizon=${HORIZON} \
  flip=${FLIP} gravity=${G} cfl=${CFL} frame_rate=${RATE} \
  solver_iterations=${SOLVER_ITERS} use_geometric=${GEOMETRIC} \
  output=false init_file=~/init_file"

echo "Parameters are:"
echo $params

./build/src/canary_launcher \
  --launch_application=${HOME}/research/fluid-simulations/projects/flip-current-lb/build/water_app.so \
  --launch_placement_algorithm=${PLACEMENT} --launch_num_worker=${NUM_WORKERS} \
  total_workers=${NUM_WORKERS} \
  init=${INIT} frames=${FRAMES} ngrid=${NGRID} \
  partition_x=${X} partition_y=${Y} partition_z=${Z} \
  ax=${AX} ay=${AY} az=${AZ} \
  processors=${P} horizon=${HORIZON} \
  flip=${FLIP} gravity=${G} cfl=${CFL} frame_rate=${RATE} \
  solver_iterations=${SOLVER_ITERS} use_geometric=${GEOMETRIC} \
  output=false  init_file=~/init_file\
  &> /dev/null
