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
      grid_solid_str_,
      particles_str_
    };
    GetReadWriteAccess(read, write);

    app_->Transform([=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << name_ << ": In initialize transform for ...";

      FlipSimulationT *simulation = SetVariables(task_context, read, write);

      // Scale gravity.
      T scaled_gravity = config_.gravity * (T(params_.global_dims[1])/128.0);
      LOG(INFO) << name_ << ": Setting gravity to " << scaled_gravity;

      // Initialize metadata and variables.
      int rank = task_context->GetPartitionId();
      simulation->InitializeMetadata(
        params_.global_dims, params_.partitions, rank,
        T(1.0), params_.nparticles, scaled_gravity, params_.flip_factor,
        params_.output_dir, params_.debug);
      simulation->InitializeVariables();
      VLOG(1) << name_ << ": Variables and metadata initialized ...";

      // Set forces. Scale forces by gravity.
      const common::Vec3<T> forces = config_.force * scaled_gravity;
      simulation->SetForces(forces);
      LOG(INFO) << name_ << ": Initialized external force "
                << common::ToString(forces);

      // Initialize sources, followed by solids and in the end fluids.
      // It is important to intiialize solids before fluids so that there is
      // no fluid where there are solids.
      // Solids includes boundaries.

      // Initialize sources.
      for (const auto &source : config_.sources) {
        Source<T> src_cpy = source;
        src_cpy.ScaleBox(params_.global_dims);
        src_cpy.ScaleVelocity(scaled_gravity);
        simulation->AddWaterSource(src_cpy);
        LOG(INFO) << name_ << ": Initialized water source "
                  << src_cpy.ToString();
      }  // for source

      // Initialize boundaries.
      const int boundary = config_.boundary;
      CHECK(boundary >= 0 && boundary <= 1) << "Expected boundary parameter "
        << " to be either 0 or 1. 0 is box boundary, 1 is tank boudnary.";
      if (boundary == 0) {
        simulation->SetBoxBoundary();
        LOG(INFO) << name_ << ": Initialized box boundary";
      } else {
        simulation->SetTankBoundary();
        LOG(INFO) << name_ << ": Initialized tank boundary";
      }  // if boundary

      // Initialize other solids, not boundary.
      for (const auto &solid : config_.solids) {
        CHECK(solid.shape->IsCube()) << "Only cuboid solids supported.";
        CubeT cube(static_cast<CubeT &>(*(solid.shape)));
        cube.Scale(params_.global_dims);
        const common::CoordBBox box = cube.get_coord_bbox();
        simulation->AddSolidCuboid(box);
        LOG(INFO) << name_ << ": Initialized solid with dimensions "
                  << common::ToString(box);
      }  // for solid

      // Finally initialize fluid --- first the default/background given by
      // fluid parameter, and then add additional fluids.
      const int init = config_.fluid;
      CHECK(init >= 0 && init <= 2) << "Expected fluid parameter one of "
        << "0, 1 or 2. 0 is for no default fluids, 1 for sphere drop, "
        << "2 for dam break.";
      if (init == 0) {
        LOG(INFO) << name_ << ": No default fluids";
      } else if (init == 1) {
        T height = T(params_.global_dims[1]/4);
        std::vector<common::Vec3<T>> centers;
        std::vector<T> radii;
        float cx = T(params_.global_dims[0]/2);
        float cy = T(params_.global_dims[1]*3/4);
        float cz = T(params_.global_dims[2]*5/8);
        float radius = T(params_.global_dims[0]/4)*0.8;
        centers.push_back(common::Vec3<T>(cx, cy, cz));
        radii.push_back(radius);
        simulation->InitializeMultipleWaterDrops(height, centers, radii);
        LOG(INFO) << name_ << ": Sphere drop: Sphere radius " << radius
                  << ", center " << cx << "," << cy << "," << cz
                  << ", reservoir height " << height;
      } else if (init == 2) {
        common::Vec3<T> start(0);
        common::Vec3<T> end;
        end.x() = T(params_.global_dims[0]/4)-1;
        end.y() = T(params_.global_dims[1]*3/4)-1;
        end.z() = T(params_.global_dims[2])-1;
        common::Cube<T> cube(start, end);
        simulation->InitializeWaterCuboid(cube, common::Vec3<T>(0));
        LOG(INFO) << name_ << ": Initialized tank with dimensions: "
                  << cube.ToString();
      }  // if init ...
      // Additional fluids.
      for (const auto &fluid : config_.additional_fluids) {
        if (fluid.shape->IsCube()) {
          CubeT cube(static_cast<CubeT &>(*(fluid.shape)));
          cube.Scale(params_.global_dims);
          cube.set_end(cube.get_end() - common::Vec3<T>(1));
          simulation->InitializeWaterCuboid(cube, fluid.velocity*scaled_gravity);
          LOG(INFO) << name_ << ": Initialized water cube with dimensions "
                    << cube.ToString();
        } else {
          CHECK(fluid.shape->IsSphere()) << "Only cuboid and spherical fluid "
                                         << "shapes supported.";
          SphereT sphere(static_cast<SphereT &>(*(fluid.shape)));
          sphere.Scale(params_.global_dims);
          simulation->InitializeWaterSphere(
              sphere.get_center(), sphere.get_radius(),
              fluid.velocity*scaled_gravity);
          LOG(INFO) << name_ << ": Initialized water sphere with center "
                    << common::ToString(sphere.get_center())
                    << " and radius " << sphere.get_radius();
        }  // if IsCube else ...
      }  // for fluid

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
  // Mark sources.
  AddSources();
  // Move particles in grid velocity.
  MoveParticles();
  SaveDebugData();  // 1
  // Delete fluid from sources/reset them for rest of the simulation.
  DeleteSources();

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
    if (*global_dt < params_.frame_step/10) {
      *global_dt = params_.frame_step/10;
      VLOG(1) << "Capping dt to " << *global_dt;
    }
    VLOG(1) << "Reduced dt = " << *global_dt;
    if (*step + *global_dt > params_.frame_step) {
      *global_dt = params_.frame_step - *step;
    }
    if (*step + *global_dt > 0.9*params_.frame_step) {
      *global_dt = params_.frame_step - *step;
    }
    assert(*global_dt > 0);
    LOG(INFO) << name_ << ": Final dt = " << *global_dt;
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

void FlipDriver::AddSources() {
  std::vector<std::string> read = {};
  std::vector<std::string> write = {
    simulation_str_, grid_marker_str_, grid_velocity_str_, particles_str_ };
  GetReadWriteAccess(read, write);
  app_->ReadAccess(*advanced_time_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": AddSources transform start ...";
    T advanced_time = task_context->ReadVariable(*advanced_time_);
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->AddContributionFromSources(advanced_time);
    ResetVariables(simulation);
    VLOG(1) << name_ << ": AddSources transform end ...";
  }, name_ + "_AddSources");
}  // AddSources

void FlipDriver::DeleteSources() {
  std::vector<std::string> read = {};
  std::vector<std::string> write = {
    simulation_str_, particles_str_ };
  GetReadWriteAccess(read, write);
  app_->ReadAccess(*advanced_time_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": DeleteSources transform start ...";
    T advanced_time = task_context->ReadVariable(*advanced_time_);
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->DeleteFluidInSources(advanced_time);
    ResetVariables(simulation);
    VLOG(1) << name_ << ": DeleteSources transform end ...";
  }, name_ + "_DeleteSources");
}  // DeleteSources

void FlipDriver::MoveParticles() {
  // Move particles (update particle positions) using grid velocity.
  std::vector<std::string> read = { grid_solid_str_, grid_velocity_str_ };
  std::vector<std::string> write = { simulation_str_, particles_str_ };
  GetReadWriteAccess(read, write);
  app_->ReadAccess(*local_dt_);
  app_->Transform([=](canary::CanaryTaskContext *task_context) {
    VLOG(1) << name_ << ": MoveParticles transform start ...";
    const T &dt = task_context->ReadVariable(*local_dt_);
    FlipSimulationT *simulation = SetVariables(task_context, read, write);
    simulation->MoveParticles(dt, true);
    // Delete particles inside solids.
    simulation->DeleteParticlesInSolids();
    ResetVariables(simulation);
    VLOG(1) << name_ << ": MoveParticles transform end ...";
  }, name_ + "_MoveParticles");

  // Exchange particles that have moved out + 1 layer ghost particles for
  // transferring particle velocities to grid correctly.
  SendParticles(2, true);
  ReceiveParticles(2);
}  // MoveParticles

void FlipDriver::ClearGridData() {
  // Clear grid data.
  std::vector<std::string> read = {};
  std::vector<std::string> write = {
    simulation_str_,
    grid_velocity_str_, grid_velocity_update_str_, grid_weight_str_,
    grid_phi_str_, grid_divergence_str_, grid_marker_str_, grid_pressure_str_ };
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
    simulation->AddGravity(dt);
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
  std::vector<std::string> read = { grid_solid_str_ };
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
