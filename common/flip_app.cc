#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>

#include "canary/canary.h"

#include "common/flip_app.h"
#include "common/definitions.h"
#include "common/utils.h"

namespace application {

/*
 * Load parameters from command line.
 */
void FlipApp::LoadParameter(const std::string &params) {
  std::stringstream ss;
  ss << params;
  {
    int ngrid;
    cereal::XMLInputArchive archive(ss);

    // Profile parameters
    LoadFlag("total_workers", profile_params_.total_workers, archive);  // REQUIRED
    LoadFlag("processors", profile_params_.processors, archive);
    LoadFlag("horizon", profile_params_.horizon, archive);
    LoadFlag("ax", profile_params_.affinity[0], archive);
    LoadFlag("ay", profile_params_.affinity[1], archive);
    LoadFlag("az", profile_params_.affinity[2], archive);

    // Simulation parameters
    LoadFlag("frames", total_frames_, archive);

    LoadFlag("partition_x", params_.main_partitions[0], archive);
    LoadFlag("partition_y", params_.main_partitions[1], archive);
    LoadFlag("partition_z", params_.main_partitions[2], archive);
    LoadFlag("ngrid", ngrid, archive);
    if (ngrid > 0) {
      params_.global_dims = common::Coord(ngrid);
    }
    LoadFlag("frame_rate", params_.frame_rate, archive);
    params_.frame_step = T(1) / T(params_.frame_rate);
    LoadFlag("cfl", params_.cfl, archive);

    LoadFlag("init", params_.init, archive);
    LoadFlag("gravity", params_.gravity, archive);
    LoadFlag("sx", params_.start[0], archive);
    LoadFlag("sy", params_.start[1], archive);
    LoadFlag("sz", params_.start[2], archive);
    LoadFlag("ex", params_.end[0], archive);
    LoadFlag("ey", params_.end[1], archive);
    LoadFlag("ez", params_.end[2], archive);
    LoadFlag("fx", params_.force[0], archive);
    LoadFlag("fy", params_.force[1], archive);
    LoadFlag("fz", params_.force[2], archive);

    LoadFlag("solver_iterations", params_.solver_max_iterations, archive);

    LoadFlag("nparticles", params_.nparticles, archive);
    LoadFlag("flip", params_.flip_factor, archive);
    LoadFlag("output_dir", params_.output_dir, archive);
    LoadFlag("output", params_.output, archive);
    LoadFlag("debug", params_.debug, archive);

    LoadFlag("use_geometric", profile_params_.use_geometric, archive);

    // When using coarse simulation, also need these parameters.
    LoadFlag("window", frames_each_time_, archive);
    LoadFlag("coarse_scale", coarse_scale_down_, archive);
    LoadFlag("coarse_partition_x", coarse_partitions_[0], archive);
    LoadFlag("coarse_partition_y", coarse_partitions_[1], archive);
    LoadFlag("coarse_partition_z", coarse_partitions_[2], archive);

    // Value of partitions defaults to main_partitions
    params_.partitions = params_.main_partitions;

    LOG(INFO) << "SIMULATION PARAMETERS: " << params;
  }
  int affinity_num = profile_params_.affinity[0] * profile_params_.affinity[1] * profile_params_.affinity[2];
  int partition_num = params_.partitions[0]*params_.partitions[1]*params_.partitions[2];
  CHECK_EQ(partition_num%affinity_num, 0);
}  // LoadParameter

/*
 * Application constructor.
 */
FlipApp::FlipApp() {
  params_.partitions = common::Coord(1);
  params_.main_partitions = common::Coord(1);
  params_.global_dims = common::Coord(64);
  params_.start = common::Vec3<T>(0);
  params_.end = common::Vec3<T>(1);
  params_.force = common::Vec3<T>(0,1,0);
  profile_params_.affinity = common::Coord(4);
  coarse_partitions_ = common::Coord(1);
  coarse_sim_ = nullptr;
  main_sim_ = nullptr;
}  // FlipApp

/*
 * Application destructor.
 */
FlipApp::~FlipApp() {}  // ~FlipApp

/*
 * Sets main simulation variables as migratable.
 */
void FlipApp::SetMigratableVariables() {
  CHECK(main_sim_ != nullptr);
  if (main_sim_->is_migratable()) {
    main_sim_->SetMigratableVariables();
  }
  if (coarse_sim_ && coarse_sim_->is_migratable()) {
    coarse_sim_->SetMigratableVariables();
  }
}  // SetMigratableVariables

/*
 * Simulation driver constructor.
 */
FlipDriver::FlipDriver(
  FlipApp *flip_app,
  const SimulationParameters &params,
  const ProfileParameters &profile_params,
  std::string name_str) : app_(flip_app), global_profile_(false),
    params_(params), profile_params_(profile_params), name_(name_str) {
  advanced_time_ = nullptr;
  global_advanced_time_ = nullptr;
  local_dt_ = nullptr;
  global_dt_ = nullptr;
  step_ = nullptr;
  simulation_ = nullptr;
  grid_velocity_ = nullptr;
  grid_velocity_update_ = nullptr;
  grid_weight_ = nullptr;
  grid_phi_ = nullptr;
  grid_divergence_ = nullptr;
  grid_pressure_ = nullptr;
  grid_marker_ = nullptr;
  grid_source_ = nullptr;
  particles_ = nullptr;
  local_distribution_ = nullptr;
  global_distribution_ = nullptr;
}  // FlipDriver

/*
 * Simulation driver destructor.
 */
FlipDriver::~FlipDriver() {
  if (advanced_time_)
    delete advanced_time_;
  if (global_advanced_time_)
    delete global_advanced_time_;
  if (local_dt_)
    delete local_dt_;
  if (global_dt_)
    delete global_dt_;
  if (step_)
    delete step_;
  if (simulation_)
    delete simulation_;
  if (grid_velocity_)
    delete grid_velocity_;
  if (grid_velocity_update_)
    delete grid_velocity_update_;
  if (grid_weight_)
    delete grid_weight_;
  if (grid_phi_)
    delete grid_phi_;
  if (grid_divergence_)
    delete grid_divergence_;
  if (grid_pressure_)
    delete grid_pressure_;
  if (grid_marker_)
    delete grid_marker_;
  if (grid_source_)
    delete grid_source_;
  if (particles_)
    delete particles_;
}  // ~FlipDriver

/*
 * Start timer.
 */
void FlipDriver::StartTimer() {
  app_->ReadAccess(*grid_phi_);
  app_->FlushBarrier([=](canary::CanaryTaskContext *task_context) {},
                         name_ + "_TimerStart");
  app_->ComputeTimerStart();
}

/*
 * Stop timer.
 */
void FlipDriver::EndTimer() {
  app_->ReadAccess(*grid_phi_);
  app_->FlushBarrier([=](canary::CanaryTaskContext *task_context) {},
                         name_ + "_TimerEnd");
  app_->ComputeTimerEnd();
}

/*
 * Declare variables.
 */
void FlipDriver::DeclareVariables() {
  LOG(INFO) << name_ << ": Setting up solver driver ...";
  const bool solver_fixed_iterations = false;
  int volume = params_.global_dims[0] * params_.global_dims[1] *
               params_.global_dims[2];
  T threshold = T(volume) * params_.solver_threshold;
  solver_ = new common::PCGDriver<ScalarGridT>(
    app_, name_,
    params_.solver_max_iterations, threshold, solver_fixed_iterations,
    common::IC, params_.partitions, params_.global_dims);

  const auto num_partitions = params_.partitions[0]*params_.partitions[1]*
                              params_.partitions[2];

  LOG(INFO) << name_ << ": Declaring time variables ...";
  advanced_time_ = new VariableHandle<T>(app_->DeclareVariable<T>(num_partitions));
  global_advanced_time_ = new VariableHandle<T>(app_->DeclareVariable<T>(1));
  local_dt_ = new VariableHandle<T>(app_->DeclareVariable<T>(num_partitions));
  global_dt_ = new VariableHandle<T>(app_->DeclareVariable<T>(1));
  step_ = new VariableHandle<T>(app_->DeclareVariable<T>(1));

  LOG(INFO) << name_ << ": Declaring simulation variables ...";
  simulation_ = new VariableHandle<FlipSimulationT>(
    app_->DeclareVariable<FlipSimulationT>(num_partitions));
  grid_velocity_ = new VariableHandle<VectorGridT>(
    app_->DeclareVariable<VectorGridT>(num_partitions));
  grid_velocity_update_ = new VariableHandle<VectorGridT>(
    app_->DeclareVariable<VectorGridT>(num_partitions));
  grid_weight_ = new VariableHandle<VectorGridT>(
    app_->DeclareVariable<VectorGridT>(num_partitions));
  grid_phi_ = new VariableHandle<ScalarGridT>(
    app_->DeclareVariable<ScalarGridT>(num_partitions));
  grid_divergence_ = new VariableHandle<ScalarGridT>(
    app_->DeclareVariable<ScalarGridT>(num_partitions));
  grid_pressure_ = new VariableHandle<ScalarGridT>(
    app_->DeclareVariable<ScalarGridT>(num_partitions));
  grid_marker_ = new VariableHandle<ScalarGridInt>(
    app_->DeclareVariable<ScalarGridInt>(num_partitions));
  grid_source_ = new VariableHandle<ScalarGridBool>(
    app_->DeclareVariable<ScalarGridBool>(num_partitions));
  particles_ = new VariableHandle<ParticlesT>(
    app_->DeclareVariable<ParticlesT>(num_partitions));

  solver_->SetVariables(grid_divergence_, grid_marker_, grid_pressure_);
  solver_->DeclareVariables();

  local_distribution_ = new VariableHandle<common::ProfileSliceLocal>(
    app_->DeclareVariable<common::ProfileSliceLocal>(num_partitions));
  // Set if distribution from all partitions are to be gathered
  // for computing partition assignment.
  if (global_profile_) {
    global_distribution_ = new VariableHandle<common::FluidDistribution>(
      app_->DeclareVariable<common::FluidDistribution>(1));
  }

  LOG(INFO) << name_ << ": Declared all variables ...";
}  // DeclareVariables

/*
 * Assign worker for global variables.
 */
void FlipDriver::AssignWorkerForGlobals(int idx) {
  std::vector<int> global = {idx};
  bool success = global_advanced_time_->RecordInitialPartitionPlacement(app_, global);
  success &= global_dt_->RecordInitialPartitionPlacement(app_, global);
  assert(success);
  success &= step_->RecordInitialPartitionPlacement(app_, global);
  if (global_profile_) {
    success &= global_distribution_->RecordInitialPartitionPlacement(app_, global);
  }
  solver_->AssignWorkerForGlobals(idx);
}  // AssignWorkerForGlobals

/*
 * Assign workers for partitioned variables.
 */
void FlipDriver::AssignWorkersForPartitions(const std::vector<int> &idxs) {
  bool success = local_dt_->RecordInitialPartitionPlacement(app_, idxs);
  success &= advanced_time_->RecordInitialPartitionPlacement(app_, idxs);
  success &= simulation_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= grid_velocity_->RecordInitialPartitionPlacement(app_, idxs);
  success &= grid_velocity_update_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= grid_weight_->RecordInitialPartitionPlacement(app_, idxs);
  success &= grid_phi_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= grid_divergence_->RecordInitialPartitionPlacement(app_, idxs);
  success &= grid_pressure_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= grid_marker_->RecordInitialPartitionPlacement(app_, idxs);
  success &= grid_source_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= particles_->RecordInitialPartitionPlacement(app_, idxs);
	success &= local_distribution_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  solver_->AssignWorkersForPartitions(idxs);
}  // AssignWorkersForPartitions

/*
 * Marks all partitioned variables for this driver as migratable.
 */
void FlipDriver::SetMigratableVariables() {
  local_dt_->SetUpdatePlacement(app_, true);
  advanced_time_->SetUpdatePlacement(app_, true);
  simulation_->SetUpdatePlacement(app_, true);
  grid_velocity_->SetUpdatePlacement(app_, true);
  grid_velocity_update_->SetUpdatePlacement(app_, true);
  grid_weight_->SetUpdatePlacement(app_, true);
  grid_phi_->SetUpdatePlacement(app_, true);
  grid_divergence_->SetUpdatePlacement(app_, true);
  grid_pressure_->SetUpdatePlacement(app_, true);
  grid_marker_->SetUpdatePlacement(app_, true);
  grid_source_->SetUpdatePlacement(app_, true);
  particles_->SetUpdatePlacement(app_, true);
  local_distribution_->SetUpdatePlacement(app_, true);
  solver_->SetMigratableVariables();
}  // SetMigratableVariables

/*
 * Initialize global distribution.
 */
void FlipDriver::InitializeGlobalProfiling() {
  CHECK(global_profile_) << "Global profiling not on for " << name_;
  // Then, global profiling stuff.
  app_->WriteAccess(*global_distribution_);;
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << " : Initializing profiling global ...";
    common::FluidDistribution *distribution = task_context->WriteVariable(
      *global_distribution_);
    distribution->Initialize(params_.main_partitions,
                             profile_params_.total_workers,
                             profile_params_.processors,
                             profile_params_.current_load_factor);
    VLOG(1) << name_ << " : Completed initializing profiling global";
  }, name_ + "InitProfilingGlobal");  // Transform
}

/*
 * Profile simulation distribution locally.
 */
void FlipDriver::ProfileDistributionLocal(bool print_rank) {
  // Profile each locally -- count number of fluid cells across all partitions.
  std::vector<std::string> read;
  std::vector<std::string> write = {
    simulation_str_, local_distribution_str_, grid_marker_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << "Profiling simulation ... ";
    double compute_time = task_context->GetComputeTime();
    double scatter_gather_time = task_context->GetScatterGatherTime();
    task_context->SetComputeTime(0);
    task_context->SetScatterGatherTime(0);
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->Profile(compute_time, scatter_gather_time,
                        profile_params_.use_geometric);
    simulation->PrintProfilingData(name_, print_rank);
    VLOG(1) << "Profiling simulation done";
  }, name_ + "_ProfileDistribution");
}  // ProfileDistributionLocal

/*
 * Gather profile globally.
 */
void FlipDriver::GatherProfileData() {
  CHECK(global_profile_) << "Global profiling not on for " << name_;
  // Send (scatter) profiled data to partition 0.
  app_->ReadAccess(*local_distribution_);
  app_->Scatter([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << "Sending profiled data to partition 0 ...";
    const common::ProfileSliceLocal &distribution =
      task_context->ReadVariable(*local_distribution_);
    struct evbuffer *buffer = evbuffer_new();
    canary::CanaryOutputArchive archive(buffer);
    archive(distribution);
    task_context->Scatter(0, canary::RawEvbuffer(buffer));
    VLOG(1) << "Sent profiled data to partition 0";
  }, name_ + "_ScatterLocalDistribution");

  // Gather profiled data and write to fluid distribution.
  app_->WriteAccess(*global_distribution_);
  app_->ReadAccess(*global_advanced_time_);
  app_->Gather([=](canary::CanaryTaskContext *task_context) -> int {
    VLOG(1) << "Receiving profiled data at partition 0 ...";
    common::FluidDistribution *distribution =
      task_context->WriteVariable(*global_distribution_);
    int num_partitions = params_.partitions[0] * params_.partitions[1] *
                         params_.partitions[2];
    EXPECT_GATHER_SIZE(num_partitions);
    std::vector<common::ProfileSliceLocal> slices(num_partitions);
    auto recv_buffer = task_context->Gather<canary::RawEvbuffer>(); 
    VLOG(1) << "Receive buffer size = " << recv_buffer.size()
              << ", expected = " << num_partitions;
    for (size_t i = 0; i < recv_buffer.size(); ++i) {
      auto &raw_buffer = recv_buffer[i];
      canary::CanaryInputArchive archive(raw_buffer.buffer);
      archive(slices[i]); 
      evbuffer_free(raw_buffer.buffer); 
    }  // for i
    distribution->AddSlice(task_context->ReadVariable(*global_advanced_time_),
                           slices);
    distribution->Print();
    VLOG(1) << "Received profiled data at partition 0";
    return 0;
  }, name_ + "_GatherDistribution");
}  // GatherProfileData

/*
 * Update palcement -- synchronize with controller.
 */
void FlipDriver::UpdatePlacement() {
  CHECK(migratable_) << "Migratable not on for " << name_;
  app_->WriteAccess(*simulation_);
  app_->ReadAccess(*local_dt_);
  app_->ReadAccess(*advanced_time_);
  app_->ReadAccess(*simulation_);
  app_->ReadAccess(*grid_velocity_);
  app_->UpdatePlacement([=](canary::CanaryTaskContext *task_context)  {
    float advanced_time = task_context->ReadVariable(*advanced_time_);
    LOG(INFO) << "Request advance time " << advanced_time <<
                 " from partition " <<
                 task_context->GetPartitionId();
    task_context->UpdatePlacement(advanced_time);
  }, name_ + "_UpdatePlacement");
}  // UpdatePlacement

/*
 * Get read/write access to supplied variables.
 */
void FlipDriver::GetReadWriteAccess(
  const std::vector<std::string> &read, const std::vector<std::string> &write) {
  for (const auto variable : read) {
    if (variable == simulation_str_)
      app_->ReadAccess(*simulation_);
    else if (variable == grid_velocity_str_)
      app_->ReadAccess(*grid_velocity_);
    else if (variable == grid_velocity_update_str_)
      app_->ReadAccess(*grid_velocity_update_);
    else if (variable == grid_weight_str_)
      app_->ReadAccess(*grid_weight_);
    else if (variable == grid_phi_str_)
      app_->ReadAccess(*grid_phi_);
    else if (variable == grid_divergence_str_)
      app_->ReadAccess(*grid_divergence_);
    else if (variable == grid_pressure_str_)
      app_->ReadAccess(*grid_pressure_);
    else if (variable == grid_marker_str_)
      app_->ReadAccess(*grid_marker_);
    else if (variable == grid_source_str_)
      app_->ReadAccess(*grid_source_);
    else if (variable == particles_str_)
      app_->ReadAccess(*particles_);
    else if (variable == local_distribution_str_)
      app_->ReadAccess(*local_distribution_);
  }
  for (const auto variable : write) {
    if (variable == simulation_str_)
      app_->WriteAccess(*simulation_);
    else if (variable == grid_velocity_str_)
      app_->WriteAccess(*grid_velocity_);
    else if (variable == grid_velocity_update_str_)
      app_->WriteAccess(*grid_velocity_update_);
    else if (variable == grid_weight_str_)
      app_->WriteAccess(*grid_weight_);
    else if (variable == grid_phi_str_)
      app_->WriteAccess(*grid_phi_);
    else if (variable == grid_divergence_str_)
      app_->WriteAccess(*grid_divergence_);
    else if (variable == grid_pressure_str_)
      app_->WriteAccess(*grid_pressure_);
    else if (variable == grid_marker_str_)
      app_->WriteAccess(*grid_marker_);
    else if (variable == grid_source_str_)
      app_->WriteAccess(*grid_source_);
    else if (variable == particles_str_)
      app_->WriteAccess(*particles_);
    else if (variable == local_distribution_str_)
      app_->WriteAccess(*local_distribution_);
  }
}  // GetReadWriteAccess

/*
 * Set read/write variables in simulation driver.
 */
FlipDriver::FlipSimulationT *FlipDriver::SetVariables(
  canary::CanaryTaskContext *task_context,
  const std::vector<std::string> &read, const std::vector<std::string> &write) {
  // Right now, just get write access over simulation variable that contains
  // simulation state.
  // It *might* be a good idea to decouple state from simulation class,
  // and make state a variable instead of simulation,
  // instantiate simulation class every time it is needed within transforms,
  // and get only read access to state when state is not updated.
  // All members of simulation that are not variables could be moved to a
  // class FlipSimulationState.
  FlipSimulationT *simulation = task_context->WriteVariable(*simulation_);
  assert(write[0] == simulation_str_);
  for (const auto variable : read) {
    if (variable == grid_velocity_str_) {
      VectorGridT &grid_velocity =
        const_cast<VectorGridT&>(task_context->ReadVariable(*grid_velocity_));
      simulation->set_grid_velocity(&grid_velocity);
    } else if (variable == grid_velocity_update_str_) {
      VectorGridT &grid_velocity_update =
        const_cast<VectorGridT&>(task_context->ReadVariable(*grid_velocity_update_));
      simulation->set_grid_velocity_update(&grid_velocity_update);
    } else if (variable == grid_weight_str_) {
      VectorGridT &grid_weight =
        const_cast<VectorGridT&>(task_context->ReadVariable(*grid_weight_));
      simulation->set_grid_weight(&grid_weight);
    } else if (variable == grid_phi_str_) {
      ScalarGridT &grid_phi =
        const_cast<ScalarGridT&>(task_context->ReadVariable(*grid_phi_));
      simulation->set_grid_phi(&grid_phi);
    } else if (variable == grid_divergence_str_) {
      ScalarGridT &grid_divergence =
        const_cast<ScalarGridT&>(task_context->ReadVariable(*grid_divergence_));
      simulation->set_grid_divergence(&grid_divergence);
    } else if (variable == grid_pressure_str_) {
      ScalarGridT &grid_pressure =
        const_cast<ScalarGridT&>(task_context->ReadVariable(*grid_pressure_));
      simulation->set_grid_pressure(&grid_pressure);
    } else if (variable == grid_marker_str_) {
      ScalarGridInt &grid_marker =
        const_cast<ScalarGridInt&>(task_context->ReadVariable(*grid_marker_));
      simulation->set_grid_marker(&grid_marker);
    } else if (variable == grid_source_str_) {
      ScalarGridBool &grid_source =
        const_cast<ScalarGridBool&>(task_context->ReadVariable(*grid_source_));
      simulation->set_grid_source(&grid_source);
    } else if (variable == particles_str_) {
      ParticlesT &particles =
        const_cast<ParticlesT&>(task_context->ReadVariable(*particles_));
      simulation->set_particles(&particles);
    } else if (variable == local_distribution_str_) {
      common::ProfileSliceLocal &profile_slice_local =
        const_cast<common::ProfileSliceLocal&>(task_context->ReadVariable(*local_distribution_));
      simulation->set_profile(&profile_slice_local);
    }
  }  // for variable : write
  for (const auto variable : write) {
    if (variable == grid_velocity_str_) {
      VectorGridT *grid_velocity =
        task_context->WriteVariable(*grid_velocity_);
      simulation->set_grid_velocity(grid_velocity);
    } else if (variable == grid_velocity_update_str_) {
      VectorGridT *grid_velocity_update =
        task_context->WriteVariable(*grid_velocity_update_);
      simulation->set_grid_velocity_update(grid_velocity_update);
    } else if (variable == grid_weight_str_) {
      VectorGridT *grid_weight =
        task_context->WriteVariable(*grid_weight_);
      simulation->set_grid_weight(grid_weight);
    } else if (variable == grid_phi_str_) {
      ScalarGridT *grid_phi =
        task_context->WriteVariable(*grid_phi_);
      simulation->set_grid_phi(grid_phi);
    } else if (variable == grid_divergence_str_) {
      ScalarGridT *grid_divergence =
        task_context->WriteVariable(*grid_divergence_);
      simulation->set_grid_divergence(grid_divergence);
    } else if (variable == grid_pressure_str_) {
      ScalarGridT *grid_pressure =
        task_context->WriteVariable(*grid_pressure_);
      simulation->set_grid_pressure(grid_pressure);
    } else if (variable == grid_marker_str_) {
      ScalarGridInt *grid_marker =
        task_context->WriteVariable(*grid_marker_);
      simulation->set_grid_marker(grid_marker);
    } else if (variable == grid_source_str_) {
      ScalarGridBool *grid_source =
        task_context->WriteVariable(*grid_source_);
      simulation->set_grid_source(grid_source);
    } else if (variable == particles_str_) {
      ParticlesT *particles =
        task_context->WriteVariable(*particles_);
      simulation->set_particles(particles);
    } else if (variable == local_distribution_str_) {
      common::ProfileSliceLocal *profile_slice_local =
        task_context->WriteVariable(*local_distribution_);
      simulation->set_profile(profile_slice_local);
    }
  }  // for variable : write
  return simulation;
}  // SetVariables

/*
 * Reset simulation variables.
 */
void FlipDriver::ResetVariables(FlipSimulationT *simulation) {
  simulation->set_grid_velocity(nullptr);
  simulation->set_grid_velocity_update(nullptr);
  simulation->set_grid_weight(nullptr);
  simulation->set_grid_phi(nullptr);
  simulation->set_grid_divergence(nullptr);
  simulation->set_grid_pressure(nullptr);
  simulation->set_grid_marker(nullptr);
  simulation->set_grid_source(nullptr);
  simulation->set_particles(nullptr);
  simulation->set_profile(nullptr);
}  // ResetVariables

const std::string FlipDriver::local_dt_str_ = "local_dt";
const std::string FlipDriver::global_dt_str_ = "global_dt";
const std::string FlipDriver::step_str_ = "step";
const std::string FlipDriver::simulation_str_ = "simulation";
const std::string FlipDriver::grid_velocity_str_ = "grid_velocity";
const std::string FlipDriver::grid_velocity_update_str_ = "grid_velocity_update";
const std::string FlipDriver::grid_weight_str_ = "grid_weight";
const std::string FlipDriver::grid_phi_str_ = "grid_phi";
const std::string FlipDriver::grid_divergence_str_ = "grid_divergence";
const std::string FlipDriver::grid_pressure_str_ = "grid_pressure";
const std::string FlipDriver::grid_marker_str_ = "grid_marker";
const std::string FlipDriver::grid_source_str_ = "grid_source";
const std::string FlipDriver::particles_str_ = "particles";
const std::string FlipDriver::local_distribution_str_ = "local_distribution";

}  // namespace application
