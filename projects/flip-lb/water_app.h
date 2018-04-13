#ifndef PROJECTS_FLIP_LB_WATER_APP_H
#define PROJECTS_FLIP_LB_WATER_APP_H

#include "common/flip_app.h"

namespace application {

class WaterApp : public FlipApp {
  public:
    virtual void Program();
    virtual void ComputeInitialPartitionPlacement(int num_workers);
};  // class WaterApp

class MainSimulationDriver : public FlipDriver {
  public:

    MainSimulationDriver(
      WaterApp *water_app, const SimulationParameters &params,
      const ProfileParameters &profile_params, std::string name = "Main")
      : FlipDriver(water_app, params, profile_params, name) {}

    // Simulation loop.
    virtual void RunSimulation(int num_frames);

};  // class MainSimulationDriver

class CoarseSimulationDriver : public FlipDriver {
  public:

    CoarseSimulationDriver(
      WaterApp *water_app, const SimulationParameters &params,
      const ProfileParameters &profile_params, std::string name = "Coarse")
      : FlipDriver(water_app, params, profile_params, name) {}

    // Simulation loop.
    virtual void RunSimulation(int num_frames);

    // Compute and send partitions to controller.
    void ComputeAndSendPartitions();

};  // class CoarseSimulationDriver

}  // namespace application

#endif  // PROJECTS_FLIP_LB_WATER_APP_H
