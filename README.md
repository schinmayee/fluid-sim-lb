# Load Balancing Fluid Simulations

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
- [OpenVDB](http://www.openvdb.org/) edited, for distributing solver and
   supporting additional types, you can download the updated version
   [here](https://github.com/schinmayee/openvdb), or use `openvdb.patch` to
   add the edits
- [Boost](https://www.boost.org/)
- [Google flags](https://gflags.github.io/gflags/)
- [Google logging library](https://github.com/google/glog)
- [Canary](https://github.com/quhang/canary), you will need to use the
  [`lb`](https://github.com/quhang/canary/tree/lb)
  branch that includes updates to use custom partition assignments computed by
  application and track simulation time at the controller.

The code has been tested with g++ 5.4, OpenVDB 4.0.1 and Boost 1.61.0 on
Ubuntu 16.04.


## Building and Running

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

3. Running:
A project will either be built as a library that can be launched over workers
using Canary, or it will be built as an executable that you can launch directly.
All generated build files will be in `build` directory. For example, the build
files for `common` will be in `common/build`, and the build files for
`projects/flip-lb` will be in `projects/flip-lb/build`.
You can view command line options for executables using `--helpshort`.
