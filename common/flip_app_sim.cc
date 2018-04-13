#include <cstdlib>
#include <iostream>
#include <string>

#include "canary/canary.h"

#include "common/flip_app.h"
#include "common/definitions.h"
#include "common/utils.h"

#define THRESH 1E-5

namespace application {

/*
 * Min of two values.
 */
template<typename T>
T min(const T value1, const T value2) {
  if (value1 < value2) {
    return value1;
  } else {
    return value2;
  }
}

/*
 * Initialize simulation.
 */
void FlipDriver::Initialize() {
  // Initialize parameters.
  app_->WriteAccess(*step_);
  app_->WriteAccess(*global_advanced_time_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": Initializing global params ...";
    T *step = task_context->WriteVariable(*step_);
    *step = 0;
    T *global_advanced_time = task_context->WriteVariable(*global_advanced_time_);
    *global_advanced_time = 0;
    VLOG(1) << name_ << ": Intiialized global params ...";
  }, name_ + "_InitGlobalParameters");
  app_->WriteAccess(*advanced_time_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": Initializing local params ...";
    T *advanced_time = task_context->WriteVariable(*advanced_time_);
    *advanced_time = 0;
    VLOG(1) << name_ << ": Intiialized local params ...";
  }, name_ + "_InitPerPartitionParameters");
  
  // Initialize simulation grid/particle variables.
  {
    // Get read/write access.
    std::vector<std::string> read;
    std::vector<std::string> write = {
      simulation_str_,
      grid_velocity_str_,
      grid_velocity_update_str_,
      grid_weight_str_,
      grid_phi_str_,
      grid_divergence_str_,
      grid_pressure_str_,
      grid_marker_str_,
      grid_source_str_,
      particles_str_
    };
    GetReadWriteAccess(read, write);

    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": In initialize transform for ...";

      FlipSimulationT *simulation = SetVariables(task_context, read, write);

      // Initialize metadata and variables.
      int rank = task_context->GetPartitionId();
      T gravity = params_.gravity;
      simulation->InitializeMetadata(
        params_.global_dims, params_.partitions, rank,
        T(1.0), params_.nparticles, params_.gravity, params_.flip_factor,
        params_.output_dir, params_.debug);
      simulation->InitializeVariables();
      VLOG(1) << name_ << ": Variables and metadata initialized ...";

      // Initialize simulation condition.
      int init = (params_.init > 10)? 1 : params_.init;
      switch(init) {
        case 1: {
          T height = T(params_.global_dims[1])/T(16);
          common::Vec3<T> center(
            T(params_.global_dims[0])/T(2), T(6)*T(params_.global_dims[1])/T(8),
            T(params_.global_dims[2])/T(2));
          T radius = 0.66*T(params_.global_dims[1])/T(4);
          simulation->InitializeWaterDrop(height, center, radius);
          VLOG(1) << name_ << ": Initialized water sphere drop ...";
          break;
        }
        case 2: {
          T width  = T(params_.global_dims[0]*3/8);
          T height = T(params_.global_dims[1]*5/8);
          T depth = T(params_.global_dims[2]/2);
          simulation->InitializeOneWayDamBreak(width, height, depth);
          VLOG(1) << name_ << ": Initialized one way dam break ...";
          break;
        }
        case 3: {
          T height = T(params_.global_dims[1])/8;
          simulation->InitializeDisturbance(height);
          VLOG(1) << name_ << ": Initialized disturbance ...";
          break;
        }
        case 4: {
          T height = T(params_.global_dims[1])/T(64);
          const int num_spheres = 8;
          std::vector<common::Vec3<T>> centers;
          std::vector<T> radii;
          for (int i = 0; i < num_spheres; ++i) {
            int x = i%2;
            int y = (i/2)%2;
            int z = i/4;
            float cx = T((2*x+1)*params_.global_dims[0]/4);
            float cy = T((2*y+2)*params_.global_dims[1]/5);
            float cz = T((2*z+1)*params_.global_dims[2]/4);
            float radius = T(params_.global_dims[0]*(num_spheres+i)/(12*num_spheres));
            centers.push_back(common::Vec3<T>(cx, cy, cz));
            radii.push_back(radius);
          }
          simulation->InitializeMultipleWaterDrops(height, centers, radii);
          break;
        }
        case 5: {
          T height = T(params_.global_dims[1])*params_.end.y();
          const int num_spheres = 4;
          std::vector<common::Vec3<T>> centers;
          std::vector<T> radii;
          for (int i = 0; i < num_spheres; ++i) {
            int x = i%2;
            int z = i/2;
            float cx = T((2*x+1)*params_.global_dims[0]/4);
            float cy = T(3*params_.global_dims[1]/4);
            float cz = T((2*z+1)*params_.global_dims[2]/4);
            float radius = T(params_.global_dims[0]/8)*1.5;
            centers.push_back(common::Vec3<T>(cx, cy, cz));
            radii.push_back(radius);
          }
          simulation->InitializeMultipleWaterDrops(height, centers, radii);
          break;
        }
        case 6: {
          T ax = T(1*params_.global_dims[0]/8);
          T ay = T(params_.global_dims[1]*6.5/8);
          T az = T(2*params_.global_dims[2]/8);
          T fx = params_.gravity;
          T fy = params_.gravity*0.1;
          T fz = params_.gravity;
          simulation->InitializeSloshingTank(0, ax, 0, ay, 0, az, fx, fy, fz);
          break;
        }
        case 7: {
          T height = T(params_.global_dims[1]);
          T width  = T(params_.global_dims[0])/2.0;
          T depth  = T(params_.global_dims[0])/2.0;
          simulation->InitializeOneWayDamBreak(width, height, depth);
          break;
        }
        case 8: {
          T height = T(params_.global_dims[1])*params_.end.y();
          std::vector<common::Vec3<T>> centers;
          std::vector<T> radii;
          float cx = T(params_.global_dims[0]/4);
          float cy = T(params_.global_dims[1]*3/4);
          float cz = T(params_.global_dims[2]*5/8);
          float radius = T(params_.global_dims[0]/4)*0.8*params_.end.x();
          centers.push_back(common::Vec3<T>(cx, cy, cz));
          radii.push_back(radius);
          simulation->InitializeMultipleWaterDrops(height, centers, radii);
          break;
        }
        case 9: {
          T ex = params_.end.x() * T(params_.global_dims[0]);
          T ey = params_.end.y() * T(params_.global_dims[1]);
          T ez = params_.end.z() * T(params_.global_dims[2]);
          LOG(INFO) << "Dam params: " << "Width " << ex << ", Height " << ey
                    << ", Depth " << ez;
          simulation->InitializeOneWayDamBreak(ex, ey, ez);
          break;
        }
        case 10: {
          T sx = params_.start.x() * T(params_.global_dims[0]);
          T sy = params_.start.y() * T(params_.global_dims[1]);
          T sz = params_.start.z() * T(params_.global_dims[2]);
          T ex = params_.end.x() * T(params_.global_dims[0]);
          T ey = params_.end.y() * T(params_.global_dims[1]);
          T ez = params_.end.z() * T(params_.global_dims[2]);
          T fx = params_.force.x() * params_.gravity;
          T fy = params_.force.y() * params_.gravity;
          T fz = params_.force.z() * params_.gravity;
          LOG(INFO) << "Tank params: Width " << sx << " to " << ex
                    << ", Height " << sy << " to " << ey
                    << ", Depth" << sz << " to " << ez
                    << "Force " << fx << "," << fy << "," << fz;
          simulation->InitializeSloshingTank(sx, ex, sy, ey, sz, ez, fx, fy, fz);
          break;
        }
        default: {
          std::cerr << name_ << ": Initialize simulation in invalid branch!!";
          VLOG(1) << "Init = " << init;
          assert(false);
        }
      };  // switch

      ResetVariables(simulation);
    }, name_ + "_Init");  // Transform
  }

  // Initialize local profiling stuff.
  {
    std::vector<std::string> read;
    std::vector<std::string> write = { simulation_str_, local_distribution_str_ };
    GetReadWriteAccess(read, write);
    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << " : Initializing profiling local ...";
      FlipSimulationT *simulation = SetVariables(task_context, read, write);
      simulation->InitializeProfiling(params_.main_partitions);
      VLOG(1) << name_ << " : Completed initializing profiling local";
    }, name_ + "InitProflingLocal");  // Transform
  }

  // Save data after initializing.
  SaveData();
  SaveDebugData();
}  // Initialize

/*
 * Run one time step.
 */
void FlipDriver::RunOneStep() {
  // Move particles in grid velocity.
  MoveParticles();
  SaveDebugData();  // 1

  // Clear grid data.
  ClearGridData();
  SaveDebugData();  // 2

  // Transfer particle velocities to grid, and save grid velocity.
  TransferParticlesToGrid();
  SaveDebugData(); // 3

  // Reinitialize phi outside.
  ReinitializePhi();
  SaveDebugData(); // 4

  SweepPhi();
  SaveDebugData(); // 5

  // Add gravity, advance grid velocity.
  AdvanceGridVelocity();
  SaveDebugData(); // 6

  // Apply boundary conditions.
  ApplyBoundaryConditions();
  SaveDebugData();  // 7
  // Make incompressible.
  MakeIncompressible();
  SaveDebugData();  // 8

  // Get velocity update.
  ComputeVelocityUpdate();
  SaveDebugData();  // 9

  // Update particle velocities using grid velocity.
  UpdateParticles();
  SaveDebugData();  // 10
}  // RunOneStep

/*
 * Save simulation data.
 */
void FlipDriver::SaveData() {
  if (!params_.output) {
    LOG(INFO) << "Will not write output data as " << name_
              << " is being run with output flag set to false.";
    return;
  }
  std::vector<std::string> read = {
    grid_velocity_str_, grid_phi_str_, particles_str_ };
  std::vector<std::string> write = { simulation_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << "Saving data ...";
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->SaveDataToFile();
    simulation->NextFrame();
  }, name_ + "_SaveData");
}  // SaveData

/*
 * Save simulation data for debugging when debug flag is on.
 */
void FlipDriver::SaveDebugData() {
  if (!params_.debug) {
    VLOG(1) << "Will not write debug data as " << name_ 
            << "is not being run with debug flag.";
    return;
  }
  std::vector<std::string> read = {
    grid_velocity_str_, grid_phi_str_, grid_pressure_str_, particles_str_ };
  std::vector<std::string> write = { simulation_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << "Saving debug data ...";
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->SaveDebugDataToFile();
  }, name_ + "_SaveDebugData");
}  // SaveDebugData

void FlipDriver::SendParticles(int ghost_width, bool send_outside_particles) {
  if (ghost_width < 0) {
    std::cerr << "Ignoring SendParticles call with ghost width " <<
                 ghost_width;
    return;
  }  // if ghost_width < 0
  app_->ReadAccess(*particles_);
  app_->Scatter([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": SendParticles start ...";
    // Get access to particles.
    const ParticlesT &particles = task_context->ReadVariable(*particles_);
    particles.SendGhostData(task_context, ghost_width, send_outside_particles);
    VLOG(1) << name_ << ": SendParticles end ...";
  });
}  // SendParticles

void FlipDriver::ReceiveParticles(int ghost_width) {
  app_->WriteAccess(*particles_);
  app_->Gather([=](canary::CanaryTaskContext *task_context) -> int {
    VLOG(1)  << name_ << ": ReceiveParticles start ...";
    // Get access to particles.
    ParticlesT *particles = task_context->WriteVariable(*particles_);
    int ret = particles->ReceiveGhostData(task_context, ghost_width);
    VLOG(1)  << name_ << ": ReceiveParticles end ...";
    return ret;
  });
}  // ReceiveParticles

template<typename GridT>
void FlipDriver::SendGhostData(
  VariableHandle<GridT> *variable, int ghost_width) {
  if (ghost_width <= 0) {
    std::cerr << "Ignoring SendGhostData call with ghost width " <<
                 ghost_width;
    return;
  }  // if ghost_width <= 0
  app_->ReadAccess(*variable);
  app_->Scatter([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": SendGhostData start ... ";
    // Get access to data.
    const GridT &data = task_context->ReadVariable(*variable);
    data.SendGhostData(task_context, ghost_width);
    VLOG(1) << name_ << ": SendGhostData end ...";
  });
}  // SendGhostData

template<typename GridT>
void FlipDriver::ReceiveGhostData(
  VariableHandle<GridT> *variable, int ghost_width) {
  app_->WriteAccess(*variable);
  app_->Gather([=](canary::CanaryTaskContext *task_context) -> int {
    VLOG(1) << name_ << ": ReceiveGhostData start ... ";
    // Get access to data.
    GridT *data = task_context->WriteVariable(*variable);
    int ret = data->ReceiveGhostData(task_context, ghost_width);
    VLOG(1) << name_ << ": ReceiveGhostData end ...";
    return ret;
  });
}  // ReceiveGhost

void FlipDriver::WhileAndComputeDt() {
  // While current time < time to advance to, execute more loops.
  app_->ReadAccess(*step_);
  app_->While([=](canary::CanaryTaskContext *task_context) -> bool {
    // Check current time in step.
    const T &step = task_context->ReadVariable(*step_);
    return (step + THRESH < params_.frame_step);
  }, name_ + "_WhileDt");

  // Compute dt.
  std::vector<std::string> read = { grid_velocity_str_ };
  std::vector<std::string> write = { simulation_str_ };
  GetReadWriteAccess(read, write);
  app_->WriteAccess(*local_dt_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    T *local_dt = task_context->WriteVariable(*local_dt_);
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    *local_dt = T(params_.cfl) * simulation->CFL();
    LOG(INFO) << name_ << ": Locally computed dt at partition " <<
                 task_context->GetPartitionId() << " = " << *local_dt;
    ResetVariables(simulation);
  }, name_ + "_LocalDt");

  // Reduce dt.
  // First scatter.
  app_->ReadAccess(*local_dt_);
  app_->Scatter([=](canary::CanaryTaskContext *task_context) {
    task_context->Scatter(0, task_context->ReadVariable(*local_dt_));
  });
  // Then gather and update current step.
  app_->WriteAccess(*global_dt_);
  app_->WriteAccess(*step_);
  app_->Gather([=](canary::CanaryTaskContext *task_context) -> int {
    int num_gathers = task_context->GetScatterParallelism();
    EXPECT_GATHER_SIZE(num_gathers);
    VLOG(1) << "Number of gathers = " << num_gathers;
    T *global_dt = task_context->WriteVariable(*global_dt_);
    T *step = task_context->WriteVariable(*step_);
    *global_dt = 1.0/params_.frame_step;
    *global_dt = task_context->Reduce(*global_dt, min<T>);
    if (*global_dt < params_.frame_step/100) {
      *global_dt = params_.frame_step/100;
      VLOG(1) << "Capping dt to " << *global_dt;
    }
    VLOG(1) << "Reduced dt = " << *global_dt;
    //*global_dt = .5;
    if (*step + *global_dt > params_.frame_step) {
      *global_dt = params_.frame_step - *step;
    }
    assert(*global_dt > 0);
    VLOG(1) << name_ << ": Final dt = " << *global_dt;
    (*step) += (*global_dt);
    VLOG(1) << name_ << ": Step = " << *step;
    return 0;
  }, name_ + "_ReduceDt");

  // Finally broadcast dt.
  app_->ReadAccess(*global_dt_);
  app_->Scatter([=](canary::CanaryTaskContext *task_context) {
    task_context->Broadcast(task_context->ReadVariable(*global_dt_));
  }, name_ + "_BroadcastDt");

  app_->TrackNeeded();
  app_->WriteAccess(*local_dt_);
  app_->Gather([=](canary::CanaryTaskContext *task_context) -> int {
    EXPECT_GATHER_SIZE(1);
    task_context->GatherSingle(task_context->WriteVariable(*local_dt_));
    return 0;
  }, name_ + "_ReceiveGlobalDt");
}  // WhileAndComputeDt

void FlipDriver::ResetDt() {
  app_->WriteAccess(*step_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    T *step = task_context->WriteVariable(*step_);
    *step = 0;
  }, name_ + "_ResetDt");
}  // ResetDt

void FlipDriver::AdvanceTime() {
  app_->ReadAccess(*global_dt_);
  app_->WriteAccess(*global_advanced_time_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << "Updating global advanced time ...";
    T *global_advanced_time = task_context->WriteVariable(*global_advanced_time_);
    *global_advanced_time += task_context->ReadVariable(*global_dt_);
    VLOG(1) << "Updated global advanced time";
  }, name_ + "_AdvanceTime");

  app_->ReadAccess(*global_advanced_time_);
  app_->Scatter([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << "Broadcasting advanced time ...";
    task_context->Broadcast(task_context->ReadVariable(*global_advanced_time_));
  }, name_ + "_BroadcastAdvancedTime");

  app_->WriteAccess(*advanced_time_);
  app_->Gather([=](canary::CanaryTaskContext *task_context) -> int {                
    VLOG(1) << "Updating local advance time ...";
    EXPECT_GATHER_SIZE(1);                                                         
    task_context->GatherSingle(task_context->WriteVariable(*advanced_time_));         
    VLOG(1) << "Updated local advance time";
    return 0;                                                                      
  }, name_ + "_UpdateLocalAdvancedTime"); 
}  // AdvanceTime

void FlipDriver::MoveParticles() {
  // Move particles (update particle positions) using grid velocity.
  std::vector<std::string> read = {grid_velocity_str_};
  std::vector<std::string> write = { simulation_str_, particles_str_ };
  GetReadWriteAccess(read, write);
  app_->ReadAccess(*local_dt_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": MoveParticles transform start ...";
    const T &dt = task_context->ReadVariable(*local_dt_);
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->MoveParticles(dt, true);
    ResetVariables(simulation);
    VLOG(1) << name_ << ": MoveParticles transform end ...";
  }, name_ + "_MoveParticles");

  // Exchange particles that have moved out + 1 layer ghost particles for
  // transferring particle velocities to grid correctly.
  SendParticles(1, true);
  ReceiveParticles(1);
}  // MoveParticles

void FlipDriver::ClearGridData() {
  // Clear grid data.
  std::vector<std::string> read = {};
  std::vector<std::string> write = {
    simulation_str_,
    grid_velocity_str_, grid_velocity_update_str_, grid_weight_str_,
    grid_phi_str_, grid_divergence_str_, grid_marker_str_, grid_source_str_,
    grid_pressure_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": ClearGridData start ...";
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->ClearGridData();
    ResetVariables(simulation);
    VLOG(1) << name_ << ": ClearGridData end ...";
  }, name_ + "_ClearGridData");
}  // ClearGridData

void FlipDriver::TransferParticlesToGrid() {
  // Transfer particle velocities to grid.
  {
    std::vector<std::string> read = { particles_str_ };
    std::vector<std::string> write = {
      simulation_str_, grid_marker_str_, grid_velocity_str_, grid_weight_str_,
      grid_phi_str_ };
    GetReadWriteAccess(read, write);
    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": TransferParticles start ...";
      FlipSimulationT *simulation = SetVariables(task_context, read, write);
      simulation->TransferParticlesToGrid();
      ResetVariables(simulation);
      VLOG(1) << name_ << ": TransferParticles end ...";
    }, name_ + "_TransferParticles");
  }

  // Apply boundary conditions.
  ApplyBoundaryConditions();

  // Save grid velocities.
  {
    std::vector<std::string> read = { grid_velocity_str_ };
    std::vector<std::string> write = {
      simulation_str_, grid_velocity_update_str_ };
    GetReadWriteAccess(read, write);
    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": SaveGridVelocities start ...";
      FlipSimulationT *simulation = SetVariables(task_context, read, write);
      simulation->SaveGridVelocities();
      ResetVariables(simulation);
      VLOG(1) << name_ << ": SaveGridVelocities end...";
    }, name_ + "_SaveGridVelocities");
  }

  // Delete outside particles.
  {
    std::vector<std::string> read = {};
    std::vector<std::string> write = {
      simulation_str_, particles_str_ };
    GetReadWriteAccess(read, write);
    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": DeleteOutsideParticles start ...";
      FlipSimulationT *simulation = SetVariables(task_context, read, write);
      simulation->DeleteOutsideParticles();
      ResetVariables(simulation);
      VLOG(1) << name_ << ": DeleteOutsideParticles end...";
    }, name_ + "_DeleteOutsideParticles");
  }
}  // TransferParticlesToGrid

void FlipDriver::ReinitializePhi() {
  std::vector<std::string> read = { grid_marker_str_ };
  std::vector<std::string> write = { simulation_str_, grid_phi_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": ReinitializePhi start ...";
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->ReinitializePhiSimple();
    ResetVariables(simulation);
    VLOG(1) << name_ << ": ReinitializePhi end ...";
  }, name_ + "_ReinitializePhi");
}  // ReinitializePhi

void FlipDriver::SweepPhi() {
  const int num_sweeps = 0;
  for (int s = 0; s < num_sweeps; ++s) {
    SendGhostData(grid_phi_, params_.cfl+1);
    ReceiveGhostData(grid_phi_, params_.cfl+1);
    std::vector<std::string> read = { grid_marker_str_ };
    std::vector<std::string> write = { simulation_str_, grid_phi_str_ };
    GetReadWriteAccess(read, write);
    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": SweepPhi start ...";
      FlipSimulationT *simulation = SetVariables(task_context, read, write);
      simulation->SweepPhi();
      ResetVariables(simulation);
      VLOG(1) << name_ << ": SweepPhi end ...";
    }, name_ + "_SweepPhi");
  }  // for s
  SendGhostData(grid_phi_, params_.cfl+1);
  ReceiveGhostData(grid_phi_, params_.cfl+1);
}  // SweepPhi

void FlipDriver::AdvanceGridVelocity() {
  // Add gravity.
  std::vector<std::string> read = {};
  std::vector<std::string> write = {
    simulation_str_, grid_velocity_str_ };
  GetReadWriteAccess(read, write);
  app_->ReadAccess(*local_dt_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": AdvanceGridVelocity start ...";
    const T &dt = task_context->ReadVariable(*local_dt_);
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->AddGravity(dt, true);
    ResetVariables(simulation);
    VLOG(1) << name_ << ": AdvanceGridVelocity end ...";
  }, name_ + "_AdvanceGridVelocity");

  // Exchange CFL layers of velocity
  SendGhostData(grid_velocity_, params_.cfl);
  ReceiveGhostData(grid_velocity_, params_.cfl);
}  // AdvanceGridVelocity

void FlipDriver::MakeIncompressible() {
  {
    std::vector<std::string> read = { grid_marker_str_, grid_velocity_str_ };
    std::vector<std::string> write = {
      simulation_str_, grid_divergence_str_ };
    GetReadWriteAccess(read, write);
    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": MakeIncompressible start ...";

      VLOG(1) << name_ << ": Computing divergence ...";
      FlipSimulationT *simulation = SetVariables(task_context, read, write);
      simulation->ComputeDivergence();
      ResetVariables(simulation);
    }, name_ + "_ComputeDivergence");
  }

  app_->TrackNeeded();
  solver_->Solve();
  app_->TrackNeeded();

  // Exchange ghost pressure.
  SendGhostData(grid_pressure_, 1);
  ReceiveGhostData(grid_pressure_, 1);

  {
    std::vector<std::string> read = { grid_marker_str_, grid_pressure_str_ };
    std::vector<std::string> write = {
      simulation_str_, grid_velocity_str_ };
    GetReadWriteAccess(read, write);
    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": Updating velocity from pressure ...";
      FlipSimulationT *simulation = SetVariables(task_context, read, write);
      simulation->UpdateVelocityFromPressure();
      ResetVariables(simulation);
      VLOG(1) << name_ << ": MakeIncompressible end ...";
    }, name_ + "_UpdateVelocityFromPressure");
  }

  // Exchange CFL layers of velocity
  SendGhostData(grid_velocity_, params_.cfl);
  ReceiveGhostData(grid_velocity_, params_.cfl);
}  // MakeIncompressible

void FlipDriver::ApplyBoundaryConditions() {
  std::vector<std::string> read = { grid_source_str_ };
  std::vector<std::string> write = {
    simulation_str_, grid_marker_str_, grid_velocity_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": ApplyBoundaryConditions start ...";
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->ApplyBoundaryConditions();
    ResetVariables(simulation);
    VLOG(1) << name_ << ": ApplyBoundaryConditions end ...";
  }, name_ + "_ApplyBoundaryConditions");
}  // ApplyBoundaryConditions

void FlipDriver::ComputeVelocityUpdate() {
  // Compute velocity update.
  std::vector<std::string> read = { grid_velocity_str_ };
  std::vector<std::string> write = {
    simulation_str_, grid_velocity_update_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": ComputeVelocityUpdate start ...";
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->GetVelocityUpdate();
    ResetVariables(simulation);
    VLOG(1) << name_ << ": ComputeVelocityUpdate end ...";
  }, name_ + "_ComputeVelocityUpdate");

  // Exchange 1 ghost layer of velocity update (dv).
  SendGhostData(grid_velocity_update_, 1);
  ReceiveGhostData(grid_velocity_update_, 1);
}  // ComputeVelocityUpdate

void FlipDriver::UpdateParticles() {
  // Combine grid velocity and velocity update, to compute particle velocities.
  std::vector<std::string> read = {};
  std::vector<std::string> write = {
    simulation_str_, grid_velocity_str_, grid_velocity_update_str_,
    particles_str_ };
  GetReadWriteAccess(read, write);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": UpdateParticles start ...";
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->UpdateParticlesFromGrid();
    ResetVariables(simulation);
    VLOG(1) << name_ << ": UpdateParticles end partition "
              << task_context->GetPartitionId() << " ...";
  }, name_ + "_UpdateParticles");
}  // UpdateParticles

template void FlipDriver::SendGhostData<FlipDriver::VectorGridT>(
  VariableHandle<FlipDriver::VectorGridT> *, int);
template void FlipDriver::ReceiveGhostData<FlipDriver::VectorGridT>(
  VariableHandle<FlipDriver::VectorGridT> *, int);

template void FlipDriver::SendGhostData<FlipDriver::ScalarGridT>(
  VariableHandle<FlipDriver::ScalarGridT> *, int);
template void FlipDriver::ReceiveGhostData<FlipDriver::ScalarGridT>(
  VariableHandle<FlipDriver::ScalarGridT> *, int);

}  // namespace application
