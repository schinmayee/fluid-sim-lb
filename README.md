# Load Balancing Fluid Simulations

Updates and fixes to this code are avaiable on Github [here](https://github.com/schinmayee/fluid-sim-lb).

The aim of this project is to load balance fluid simulations over dynamic data
structures. This project includes code to distribute and load balance
simulations using different strategies, and an implementation for a basic
[FLIP](https://dl.acm.org/citation.cfm?id=1073298) simulations that uses
[OpenVDB](http://www.openvdb.org/) for storing grid fields and particles.


## Dependencies

You will need the following installed and running, in order to build and run
this project:
- g++
- GNU Make
- [Google flags](https://gflags.github.io/gflags/) needed only for code that
  does particle to signed distance conversion
- [Google logging library](https://github.com/google/glog) needed only for code
  that does particle to signed distance conversion
- [OpenVDB](http://www.openvdb.org/) edited, for distributing solver and
   supporting additional types, you can download the updated version
   [here](https://github.com/schinmayee/openvdb), or use `openvdb.patch` to
   add the edits
- [Boost](https://www.boost.org/)
- [Canary](https://github.com/quhang/canary), you will need to use the
  [`lb`](https://github.com/quhang/canary/tree/lb)
  branch that includes updates to use custom partition assignments computed by
  application and track simulation time at the controller.

The code has been tested with g++ 5.4, OpenVDB 4.0.1 and Boost 1.61.0 on
Ubuntu 16.04.


## Code Structure

* `common` contains FLIP simulation code and code to compute partition to
worker assignments

* `projects/flip-current-lb` contains a FLIP application/driver that uses
geometric/reactive load balancing

* `projects/flip-lb` contains a FLIP application/driver that uses
speculative load balancing

* `projects/particles-to-levelset` contains code to convert particles to
signed distance

* `scripts` contains scripts to launch experiments on Google cloud, collect
results and parse log files

* `scripts/experiments` contains sample configuration files for the Google
cloud scripts


## Building

Building and running code involves the following steps:

1. Editing dependency paths:
Start by copying `makeinclude/paths-template.mk` to `makeinclude/paths.mk`.
Edit paths for dependencies in `makeinclude/paths.mk` to point to the correct
paths on your machine.

2. Build the code:
You can build the entire code by issuing `make` in the root directory.
If you wish to build only parts of the code, you can issue `make` in the
respective directory. Note that you need to build `common` to build most
projects.

A project will either be built as a library that can be launched over workers
using Canary, or it will be built as an executable that you can launch directly.
The fluid simulations `fluid-lb` amd `fluid-current-lb` will be built as
libraries, that you can launch with Canary workers.
All generated build files will be in `build` directory. For example, the build
files for `common` will be in `common/build`, and the build files for
`projects/flip-lb` will be in `projects/flip-lb/build`.
You can view command line options for executables using `--helpshort`.


## Running

1. Launch Canary controller from Canary build directory:
```
./src/canary_controller --v=1 &
```

2. Launch Canary workers and specify the number of worker threads for each
worker, from the Canary build directory.
You can use `--helpshort` to view command line options, and launch 2 workers
each with 8 threads as follows:
```
./src/canary_worker --worker_service 40001 --worker_threads 8 &
./src/canary_worker --worker_service 40002 --worker_threads 8 &
```

3. Launch application from the Canary build directory.
For instance, you can launch `flip-current-lb` over 2 workers, with 2X2X2
partitioning over a 128-cubed grid, and run a sphere drop simulation as follows:
```
./src/canary_launcher --launch_application=${PATH_TO_THIS_PROJECT}/projects/flip-current-lb/build/water_app.so --launch_placement_algorithm=application --launch_num_worker=2 total_workers=2 partition_x=2 partition_y=2 partition_z=2 ax=2 ay=2 az=2 frames=200 ngrid=128 frame_rate=30 cfl=1 init_file=configs/tests/sphere-drop.txt
```
The above command will run sphere drop for 200 frames and generate output files
in directories `output_0` to `output_7` (since there are 8 partitions).

4. Convert particles from step 3 to signed distance for visualizing using
OpenVDB viewer, use `--helpshort` for an explanation of command line arguments:
```
mkdir sphere_drop
${PATH_TO_THIS_PROJECT}/projects/particles-to-levelset/build/convert -start 0 -end 200 -input ${PATH_TO_OUTPUT_DIRECTORIES_FROM_STEP_3}/output -output sphere_drop/levelset -partitions 8 --alsologtostderr --radius 1.6 --voxel 1.0 --filter gaussian
```

# Parameters

The following command line parameters are available for running a FLIP application:
* `total_workers`: specify number of workers available
* `partition_x`, `partition_y`, `partition_z` : number of partitions along `x`, `y` and `z` dimension
* `ngrid` : size of grid in each dimension
* `ax`, `ay`, `az` : geometric partitioning/affinity parameters
* `use_geometric` : use geometric partitioning
* `fx`, `fy`, `fz` : forces
* `gravity`: force due to gravity
* `flip` : amount of FLIP contribution (vs PIC)
* `nparticles` : number of particles per cell
* `cfl` : CFL number
* `config` : simulation configuration (see `flip_app_sim.cc`)
* `boundary` : boundary condition (see `flip_app_sim.cc`)
* `init_file` : additional fluids, sources, solids (see `configs` for examples)
* `window` : number of time steps between two partition reassignments (frequency of load balancing)
* `horizon` : number of steps to run the coarse simulation ahead by
* `coarse_scale` : (high resolution simulation `ngrid`)/(low resolution simulation `ngrid`)
* `coarse_partition_x`, `coarse_partition_y`, `coarse_partition_z` : partitioning for low resolution (speculative) simulation
