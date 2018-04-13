#include <cstdlib>
#include <iostream>
#include <openvdb/openvdb.h>
#include <string>
#include <tbb/task_scheduler_init.h>

#include "canary/canary.h"
#include "common/flip_app.h"

#include "projects/flip-current-lb/water_app.h"

namespace application {

/*
 * Canary application program to stage.
 */
void WaterApp::Program() {
  LOG(INFO) << "Running program ...";

  // Initialize tbb and openvdb.
  tbb::task_scheduler_init tbb_init(1);
  openvdb::initialize();

  main_sim_ = new SimulationDriver(this, params_, profile_params_, "Main");
  main_sim_->set_global_profile(true);
  main_sim_->set_migratable(true);

  // Declare simulation variables and metadata.
  main_sim_->DeclareVariables();
  LOG(INFO) << "DeclareVariables done";

  // Initialize main simulation.
  main_sim_->Initialize(); 
  main_sim_->InitializeGlobalProfiling(); 

  LOG(INFO) << "Will loop " << total_frames_ << " times ...";
  main_sim_->RunSimulation(total_frames_);

  LOG(INFO) << "Program staging done";
}  // Program

/*
 * Assign partitions to workers which canary expects when the simulation starts.
 * After program staging/initialization is done at controller, we will reassign
 * partitions after coarse simulation initialization is done.
 */
void WaterApp::ComputeInitialPartitionPlacement(int num_workers) {
  LOG(INFO) << "Compute initial partition placement start ...";
  // Assign workers for main simulation variables.
  const auto num_partitions = params_.partitions[0]*params_.partitions[1]*
                              params_.partitions[2];
  CHECK_EQ(num_partitions%num_workers, 0);
  const auto partitions_per_worker = num_partitions/(num_workers);
  std::vector<int> assigned_workers;
  for (int p = 0; p < num_partitions; ++p) {
    int w = (p/partitions_per_worker);
    assigned_workers.push_back(w);
  }
  main_sim_->AssignWorkerForGlobals(0);
  main_sim_->AssignWorkersForPartitions(assigned_workers);
  LOG(INFO) << "Compute initial partition placement done";
}  // ComputeInitialPartitionPlacement

/*
 * Run simulation loop.
 */
void SimulationDriver::RunSimulation(int num_frames) {
  app_->Loop(num_frames, name_+"_Loop");
  WhileAndComputeDt();

  ComputeAndSendPartitions();  // compute and send partitions every horizon steps
  UpdatePlacement();
  StartTimer();
  RunOneStep();
  EndTimer();
  ProfileDistributionLocal(/*print_rank*/ true);
  GatherProfileData();
  AdvanceTime();
  app_->EndWhile();

  ResetDt();
  SaveData();
  app_->EndLoop();
} // RunSimulation

/*
 * Compute partitioning using profile for previous step.
 */
void SimulationDriver::ComputeAndSendPartitions() {
  CHECK(global_profile_);
  app_->WriteAccess(*global_distribution_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    common::FluidDistribution *distribution =
      task_context->WriteVariable(*global_distribution_);
    // Use profile to compute partitioning.
    distribution->ComputeFixedWindowCurrentPartitionPlacement(
      0, "greedy", profile_params_.horizon, profile_params_.affinity);
    // Send partitions to controller.
    task_context->SendComputedPartitionHistory(
      distribution->GetPartitionHistory());
    // Clear profile data.
    distribution->Clear();
    // Clear partition history.
    distribution->ClearPartitionHistory();
  });
}  // ComputeAndSendPartitions

}  // namespace application

REGISTER_APPLICATION(application::WaterApp);
