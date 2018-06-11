#ifndef PROJECTS_FLIP_LB_FLIP_APP_H
#define PROJECTS_FLIP_LB_FLIP_APP_H

#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "canary/canary.h"
#include "common/distribution.h"
#include "common/flip_simulation.h"
#include "common/metagrid.h"
#include "common/particles.h"
#include "common/primitives.h"
#include "common/scalar_grid.h"
#include "common/solver_driver.h"
#include "common/sources.h"
#include "common/vector_grid.h"

namespace application {

typedef float T;
class FlipDriver;
typedef common::Shape<T> ShapeT;
typedef common::Cube<T> CubeT;
typedef common::Sphere<T> SphereT;

struct ProfileParameters {
  int total_workers = 0;  // actual number of workers, MUST PROVIDE
  int processors = 1;  // number of processors per worker
  // Horizon for load balancing.
  int horizon = 30;
  // Affinity parameters.
  common::Coord affinity; 
  T current_load_factor = 0.0;
  // Use geometric partitioning if use_geometric is set to true
  bool use_geometric = false;
};  // struct LoadBalanceParameters

struct SimulationParameters {
  common::Coord global_dims;  // number of cells in each dimension
  common::Coord partitions;  // number of partitions in each dimension
                             // equals partitions for main or coarse sim
  common::Coord main_partitions;  // number of partitions for main sim

  int frame_rate = 30;  // fps
  T frame_step = 1.0/30;  //  invserse fps (time)
  int cfl = 3;  // cfl value to use

  std::string init_file;

  int solver_max_iterations = 1500;
  T solver_threshold = 0.5*1e-8;

  int nparticles = 8;  // number of particles in each fluid cell
  T flip_factor = 0.8;  // flip factor -- pic vs flip
  std::string output_dir = "output";  // output directory
  bool output = false;  // whether to save output
  bool debug = false;  // debug mode
};  // struct SimulationParameters

struct Fluid {
  std::shared_ptr<ShapeT> shape;
  common::Vec3<T> velocity;
};  // Fluid

struct Solid {
  std::shared_ptr<ShapeT> shape;
};  // Solid

struct SimulationConfiguration {
  int fluid = 1;
  int boundary = 0;
  // Default gravity value for 32-cubed, scaled as per dimensions during
  // actual initialization of FlipSimulation.
  T gravity = 9.8;
  common::Vec3<T> force;
  std::vector<Fluid> additional_fluids;
  std::vector< Source<T> > sources;
  std::vector<Solid> solids;
};  // struct SimulationConfiguration

class FlipApp : public canary::CanaryApplication {
  friend class FlipDriver;
  public:
    // FlipApp
    FlipApp();
    ~FlipApp();

    // Program needs to be implemented by child class.
    virtual void Program() = 0;
    // Initial mapping for partitions -- also implemented by child class.
    // This runs after program, at the controller.
    virtual void ComputeInitialPartitionPlacement(int num_workers) = 0;
    // Load parameters for the program -- runs everywhere, begore program.
    virtual void LoadParameter(const std::string &params);
    // Set which variables may be migrated.
    virtual void SetMigratableVariables();

    // Parse simulation configuration
    // (sources, solids, additional fluid, forces).
    void ParseSimulationConfiguration(const std::string init_file);

  protected:
    // Simulation parameters.
    SimulationParameters params_;
    // Profile parameters.
    ProfileParameters profile_params_;
    // Simulation configuration to use.
    SimulationConfiguration config_;
    // Number of frames to simulate.
    int total_frames_ = 60;

    // Coarse simulation scale down factor.
    int coarse_scale_down_ = 8;
    common::Coord coarse_partitions_;
    // Number of frames to run coarse/main simulation at a time.
    int frames_each_time_ = 1;
    FlipDriver *coarse_sim_;  // set to nullptr when not used
    FlipDriver *main_sim_ ;
};  // class FlipApp

class FlipDriver {
  public:
    const static int N1 = 3;
    const static int N2 = 3;
    typedef common::ScalarGrid<N1, N2, T> ScalarGridT;
    typedef common::ScalarGrid<N1, N2, int32_t> ScalarGridInt;
    typedef common::ScalarGrid<N1, N2, bool> ScalarGridBool;
    typedef common::VectorGrid<N1, N2, T> VectorGridT;
    typedef common::Particles<N1, N2, T> ParticlesT;
    typedef FlipSimulation<N1, N2, T> FlipSimulationT;

    template<typename VariableType>
    using VariableHandle = canary::CanaryApplication::VariableHandle<VariableType>;

    // FlipDriver
    FlipDriver() = delete;
    FlipDriver(FlipApp *flip_app,
               const SimulationParameters &params,
               const ProfileParameters &profile_params,
               const SimulationConfiguration &config,
               std::string name_str = "unnamed");
    ~FlipDriver();

  protected:
    // FlipApp handle.
    FlipApp *app_;
    // Simulation parameters.
    SimulationParameters params_;
    // Profiling parameters.
    ProfileParameters profile_params_;
    // Simulation configuration.
    SimulationConfiguration config_;

    // Debug name
    std::string name_ = "unnamed";

    // Whether to gather profiling information, should be true if
    // using this simulation to compute partitioning.
    bool global_profile_ = false;
    // Whether the partitions of this simulation are migratable.
    bool migratable_ = false;

  public:
    // Name accessors.
    void set_name(std::string debug_name) { name_ = debug_name; }
    std::string name() { return name_; }
    // Flag accessors.
    void set_global_profile(bool profile) {
      global_profile_ = profile;
    }
    void set_migratable(bool m) {
      migratable_ = m;
    }
    bool is_migratable() const {
      return migratable_;
    }

    // Start and stop timer.
    void StartTimer();
    void EndTimer();

    // Declare simulation variables.
    void DeclareVariables();
    // Assign workers to partitions.
    void AssignWorkerForGlobals(int idx);
    void AssignWorkersForPartitions(const std::vector<int> &idxs);
    // Set partitioned variables of this simulation as migratable.
    void SetMigratableVariables();

    // Initialize global profiling.
    void InitializeGlobalProfiling();
    // Profile coarse simulation.
    void ProfileDistributionLocal(bool print_rank);
    void GatherProfileData();
    // Advance epoch -- main simulation sends request to advance epoch and
    // waits for response from controller.
    void UpdatePlacement();

    // Run post initialization simulation -- to be implemented by child class.
    virtual void RunSimulation(int num_frames) = 0;
    // Initialize simulation.
    void Initialize();
    // Run one step.
    void RunOneStep();
    // Reset dt/step values for next step.
    void ResetDt();

  protected:
    // Get read and write access.
    void GetReadWriteAccess(
      const std::vector<std::string> &read,
      const std::vector<std::string> &write);
    // Fetch simulation variable. Fetch other variables too, and set the
    // corresponding pointers in simulation.
    FlipSimulationT *SetVariables(
      canary::CanaryTaskContext *task_context,
      const std::vector<std::string> &read,
      const std::vector<std::string> &write);
    void ResetVariables(FlipSimulationT *simulation);

    // Send and receive data.
    void SendParticles(int ghost_width, bool send_outside_particles);
    void ReceiveParticles(int ghost_width);
    template<typename GridT>
    void SendGhostData(VariableHandle<GridT> *variable, int ghost_width);
    template<typename GridT>
    void ReceiveGhostData(VariableHandle<GridT> *variable, int ghost_width);

    // Save data.
    void SaveData();
    void SaveDebugData();  // does nothing if app is not in debug mode

    // Initialize different configurations.

    // While current time < time to advance to, compute dt and update time.
    void WhileAndComputeDt();
    // Advance time.
    void AdvanceTime();
    // Add source particles/mark cells and set velocity.
    void AddSources();
    // Delete fluid from sources.
    void DeleteSources();
    // Move particles in grid velocity.
    void MoveParticles();
    // Clear grid data.
    void ClearGridData();
    // Transfer particle velocities to grid and save grid velocity.
    void TransferParticlesToGrid();
    // Reinitialize phi .
    void ReinitializePhi();
    void SweepPhi();
    // Add gravity and update grid velocity.
    void AdvanceGridVelocity();
    // Enforce incompressibility by makeing fluid velocity divergence free.
    void MakeIncompressible();
    // Apply boundary conditions.
    void ApplyBoundaryConditions();
    // Compute velocity update.
    void ComputeVelocityUpdate();
    // Update particle velocities using grid velocity.
    void UpdateParticles();

  protected:
    // Handles to variables.
    VariableHandle<T> *advanced_time_ = nullptr;
    VariableHandle<T> *global_advanced_time_ = nullptr;
    VariableHandle<T> *local_dt_ = nullptr;
    VariableHandle<T> *global_dt_ = nullptr;
    VariableHandle<T> *step_ = nullptr;
    VariableHandle<FlipSimulationT> *simulation_ = nullptr;
    VariableHandle<VectorGridT> *grid_velocity_ = nullptr;
    VariableHandle<VectorGridT> *grid_velocity_update_ = nullptr;
    VariableHandle<VectorGridT> *grid_weight_ = nullptr;
    VariableHandle<ScalarGridT> *grid_phi_ = nullptr;
    VariableHandle<ScalarGridT> *grid_divergence_ = nullptr;
    VariableHandle<ScalarGridT> *grid_pressure_ = nullptr;
    VariableHandle<ScalarGridInt> *grid_marker_ = nullptr;
    VariableHandle<ScalarGridInt> *grid_solid_ = nullptr;
    VariableHandle<ScalarGridBool> *grid_source_ = nullptr;
    VariableHandle<ParticlesT> *particles_ = nullptr;
    // Profiling variables
    VariableHandle<common::ProfileSliceLocal> *local_distribution_ = nullptr;
    VariableHandle<common::FluidDistribution> *global_distribution_ = nullptr;

    // String names for variables.
    const static std::string local_dt_str_;
    const static std::string global_dt_str_;
    const static std::string step_str_;
    const static std::string simulation_str_;
    const static std::string grid_velocity_str_;
    const static std::string grid_velocity_update_str_;
    const static std::string grid_weight_str_;
    const static std::string grid_phi_str_;
    const static std::string grid_divergence_str_;
    const static std::string grid_pressure_str_;
    const static std::string grid_marker_str_;
    const static std::string grid_solid_str_;
    const static std::string grid_source_str_;
    const static std::string particles_str_;
    const static std::string local_distribution_str_;

    // Solver
    common::PCGDriver<ScalarGridT> *solver_;
};  // class FlipDriver

}  // namespace application

#endif  // PROJECTS_FLIP_LB_FLIP_APP_H
