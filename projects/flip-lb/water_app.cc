#include <cstdlib>
#include <iostream>
#include <openvdb/openvdb.h>
#include <string>
#include <tbb/task_scheduler_init.h>

#include "canary/canary.h"
#include "common/flip_app.h"

#include "projects/flip-lb/water_app.h"

namespace application {

/*
 * Canary application program to stage.
 */
void WaterApp::Program() {
  LOG(INFO) << "Running program ...";

  // Initialize tbb and openvdb.
  tbb::task_scheduler_init tbb_init(1);
  openvdb::initialize();

  // Parameters for coarse simulation.
  SimulationParameters coarse_params = params_;
  for (int d = 0; d < 3; ++d) {
    coarse_params.global_dims[d] = params_.global_dims[d]/coarse_scale_down_;
  }
  coarse_params.partitions = coarse_partitions_;
  int coarse_cfl = std::ceil(params_.cfl * 1.0 / coarse_scale_down_);
  coarse_params.cfl = coarse_cfl;
  coarse_params.gravity = params_.gravity/T(coarse_scale_down_);
  coarse_params.solver_max_iterations = params_.solver_max_iterations/coarse_scale_down_ + 1;
  coarse_params.output = false;
  coarse_params.debug = false;

  // Coarse and main simulation.
  coarse_sim_ = new CoarseSimulationDriver(
      this, coarse_params, profile_params_, "Coarse");
  coarse_sim_->set_global_profile(true);
  main_sim_ = new MainSimulationDriver(this, params_, profile_params_, "Main");
  main_sim_->set_migratable(true);

  // Declare simulation variables and metadata.
  coarse_sim_->DeclareVariables();
  main_sim_->DeclareVariables();
  LOG(INFO) << "DeclareVariables done";

  // Initialize coarse simulation and run it for window number of frames.
  coarse_sim_->Initialize(); 
  coarse_sim_->InitializeGlobalProfiling(); 
  coarse_sim_->RunSimulation(frames_each_time_ + profile_params_.horizon);

  // Now initialize main simulation and run it for window number of frames.
  main_sim_->Initialize();
  main_sim_->RunSimulation(frames_each_time_);

  // Run coarse and main simulation in a loop.
  int loop_count = total_frames_/frames_each_time_-1;
  LOG(INFO) << "Will loop " << loop_count << " X " << frames_each_time_ <<
               " times ...";
  Loop(loop_count, "Outer");
  coarse_sim_->RunSimulation(frames_each_time_);
  main_sim_->RunSimulation(frames_each_time_);
  EndLoop();

  LOG(INFO) << "Program staging done";
}  // Program

/*
 * Assign partitions to workers which canary expects when the simulation starts.
 * After program staging/initialization is done at controller, we will reassign
 * partitions after coarse simulation/before main simulation starts.
 */
void WaterApp::ComputeInitialPartitionPlacement(int num_workers) {
  LOG(INFO) << "Compute initial partition placement start ...";

  // Assign worker for coarse simulation variables.
  int coarse_num_partitions =
      coarse_partitions_[0] * coarse_partitions_[1] * coarse_partitions_[2];
  std::vector<int> coarse_workers(coarse_num_partitions, 0);
  coarse_sim_->AssignWorkerForGlobals(0);
  coarse_sim_->AssignWorkersForPartitions(coarse_workers);

  // Assign workers for main simulation variables.
  CHECK_GE(num_workers, 2) << "Expected to run with at least 2 workers";
  const auto num_partitions = params_.partitions[0]*params_.partitions[1]*
                              params_.partitions[2];
  CHECK_EQ(num_partitions%(num_workers-1), 0);
  const auto partitions_per_worker = num_partitions/(num_workers-1);
  std::vector<int> assigned_workers;
  for (int p = 0; p < num_partitions; ++p) {
    int w = (p/partitions_per_worker)+1;
    LOG(INFO) << "Adding " << w;
    assigned_workers.push_back(w);
  }
  main_sim_->AssignWorkerForGlobals(1);
  main_sim_->AssignWorkersForPartitions(assigned_workers);
  LOG(INFO) << "Compute initial partition placement done";
}  // ComputeInitialPartitionPlacement

/*
 * Main simulation loop.
 */
void MainSimulationDriver::RunSimulation(int num_frames) {
  app_->Loop(num_frames, name_+"_Loop");
  WhileAndComputeDt();

  UpdatePlacement();  // wait till controller sends proceed message.
  StartTimer();
  RunOneStep();
  EndTimer();
  ProfileDistributionLocal(/*print_rank*/ true);
  AdvanceTime();
  app_->EndWhile();

  ResetDt();
  SaveData();
  app_->EndLoop();
}  // RunMainSimulation

/*
 * Coarse simulation loop.
 */
void CoarseSimulationDriver::RunSimulation(int num_frames) {
  app_->Loop(num_frames, name_+"_Loop");
  WhileAndComputeDt();

  StartTimer();
  RunOneStep();
  EndTimer();
  ProfileDistributionLocal(false);
  GatherProfileData();
  AdvanceTime();
  ComputeAndSendPartitions();  // compute and send partitions every horizon steps
  app_->EndWhile();

  ResetDt();
  app_->EndLoop();
} // RunSimulation

/*
 * Compute partitioning using speculative profile for coarse simulation.
 */
void CoarseSimulationDriver::ComputeAndSendPartitions() {
  CHECK(global_profile_);
  app_->WriteAccess(*global_distribution_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    common::FluidDistribution *distribution =
      task_context->WriteVariable(*global_distribution_);
    // Use profile to compute partitioning.
    distribution->ComputeFixedWindowOraclePartitionPlacement(
      1, "", profile_params_.horizon, profile_params_.affinity);
    // Send partitions to controller.
    task_context->SendComputedPartitionHistory(
      distribution->GetPartitionHistory());
    distribution->ClearPartitionHistory();
  });
}  // ComputeAndSendPartitions

}  // namespace application

REGISTER_APPLICATION(application::WaterApp);
