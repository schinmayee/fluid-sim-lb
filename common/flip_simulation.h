#ifndef PROJECTS_FLIP_COMMON_WATER_SIMULATION_H
#define PROJECTS_FLIP_COMMON_WATER_SIMULATION_H

#include <cereal/archives/xml.hpp>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "canary/canary.h"
#include "common/definitions.h"
#include "common/distribution.h"
#include "common/metagrid.h"
#include "common/particles.h"
#include "common/scalar_grid.h"
#include "common/vector_grid.h"

namespace application {

/*
 * FlipSimulation class implements simulation methods.
 */
template<common::Index N1, common::Index N2, typename T = float>
class FlipSimulation {
	public:
		FlipSimulation();
    
    typedef common::VectorGrid<N1, N2, T> VectorGridT;
    typedef common::ScalarGrid<N1, N2, T> ScalarGridT;
    typedef common::ScalarGrid<N1, N2, int32_t> ScalarGridInt;
    typedef common::ScalarGrid<N1, N2, bool> ScalarGridBool;
    typedef common::Particles<N1, N2, T> ParticlesT;
    typedef typename ParticlesT::ParticleData ParticleData;

	public:
    // Don't initialize data -- only initialize parameters.
    void InitializeMetadata(
      common::Coord global_dims, common::Coord partitions,
      int rank, T voxel_len = 1, int num_particles_per_cell = 8,
      T gravity = 9.8, T flip_factor = 1, std::string output_dir = "output",
      bool debug_flag = false);

    // Initialize all data grids metadata.
    void InitializeVariables();

    // Initialize profiling data.
    void InitializeProfiling(common::Coord main_partitions);

		// These simulation methods are to be used for doing a flip simulation.

    // Initialize water drop, with existing reservoir of given height, and
    // a sphere at given center and radius.
    void InitializeWaterDrop(T height, common::Vec3<T> center, T radius);

    // Serialization and deserialization methods using archive.
    template<class Archive> void save(Archive &ar) const {
      ar(common::kMagicNumber);
      assert(profile_ == nullptr);
      ar(metagrid_);
      ar(voxel_len_, num_particles_per_cell_, gravity_);
      ar(fx_, fy_, fz_, config_);
      ar(flip_factor_, frame_, debug_, debug_step_);
      ar(phi_max_);
      ar(output_dir_, debug_output_dir_);
      ar(common::kMagicNumber);
    }  // save

    template<class Archive> void load(Archive &ar) {
      int magic_number = 0;
      ar(magic_number);
      CHECK_EQ(magic_number, common::kMagicNumber);
      ar(metagrid_);
      ar(voxel_len_, num_particles_per_cell_, gravity_);
      ar(fx_, fy_, fz_, config_);
      ar(flip_factor_, frame_, debug_, debug_step_);
      ar(phi_max_);
      ar(output_dir_, debug_output_dir_);
      magic_number = 0;
      ar(magic_number);
      CHECK_EQ(magic_number, common::kMagicNumber);
      VLOG(1) << "** Deserialized water simulation!!";
      assert(profile_ == nullptr);
      assert(abs(voxel_len_- 1.0) < 1e-2);
      SetUpOutputAndLoggingFiles();
    }  // load

    // Initialize one way dam break.
    void InitializeOneWayDamBreak(T lx, T ly, T lz);

    // Initialize disturbance.
    void InitializeDisturbance(T height);

    // Initialize reservoir.
    void InitializeReservoir(T height);

    // Initialize multiple water drops.
    void InitializeMultipleWaterDrops(
      T height, const std::vector<common::Vec3<T>> &centers,
      const std::vector<T> &radii);

    // Initialize sloshing tank.
    void InitializeSloshingTank(T mx, T ax, T my, T ay, T mz, T az, T fx, T fy, T fz);

    // Add waterfall source.
    void AddWaterFallSource(T height, T thickness);

    // Clear grid data.
    void ClearGridData();

    // Number of active cells and particles.
    inline size_t FluidCellCount() {
      return grid_marker_->data()->activeVoxelCount();
    }
    inline size_t ParticleCount() {
      return particles_->size();
    }

    // Profile -- record profiling data.
    void Profile(float compute_time, float scatter_gather_time, bool use_geometric);
    void PrintProfilingData(std::string prefix, bool print_rank);
    void ClearProfiledData();

		// CFL time step determines the largest time step that we can take
		// according to the CFL condition.
		T CFL();
		// Save current grid velocities, for computing update to grid velocities.
		void SaveGridVelocities();
    // Move particles in grid velocity.
    void MoveParticles(T dt, bool defragment);
		// Update grid velocity by adding gravity.
		void AddGravity(T dt, bool additional);
		// Advect phi.
		void AdvectPhi();

		// Recompute signed distance using fast sweeping or marching, to keep the
		// signed distance property. Phi is only useful for sweeping velocity.
		void ReinitializePhi();
    // Set distance inside fluid to be a constant = -voxel_len_/2.
		void ReinitializePhiSimple();
    // Sweep phi everywhere.
    void SweepPhi();
		// Extend velocity beyond the fluid, needed for correctly advecting phi
		// or for interpolating to particles correctly.
		void ExtendVelocity();
		// Apply boundary conditions (to velocity -- solid walls, air etc).
		void ApplyBoundaryConditions();
    // Compute divergence.
		void ComputeDivergence();
    // Update velocity using pressure.
    void UpdateVelocityFromPressure();
		// Get velocity update for grid (delta v).
		void GetVelocityUpdate();
    // Extend velocity into ghost regions.
    void ExtendGridVelocityIntoBoundaryRegion();
		// Transfer particle velocities to grid.
		void TransferParticlesToGrid();
		// Update particles using grid.
		void UpdateParticlesFromGrid();
    // Resample particles.
    void ResampleParticles();
    // Delete outside particles.
    void DeleteOutsideParticles();
    // Add contribution from sources.
    void SetContributionFromSources();

		// Write data to file.
		void SaveDataToFile();
		void SaveDebugDataToFile();
		// Advance by these many frames.
		void NextFrame();

    std::string log_file() { return log_file_; }

	private:
    common::MetaGrid metagrid_;  // partition information
		T voxel_len_;  // length of each voxel
    int num_particles_per_cell_;  // average number of particles per cell
		T gravity_;  // value for gravity
    T fx_ = 0;
    T fy_ = 0;
    T fz_ = 0;
    T flip_factor_;  // flip factor for combining grid velocity and update
		std::string output_dir_;  // output directory
		std::string debug_output_dir_;  // output directory for debug output
		bool debug_;  // debug mode
    int config_ = 0;
    int disturbance_num_ = 4;
    int applied_disturbance_ = 0;

		// Prefix filenames for saved data with these.
		std::string grid_velocity_prefix_;
		std::string grid_phi_prefix_;
		std::string particle_prefix_;
		std::string debug_grid_velocity_prefix_;
		std::string debug_grid_phi_prefix_;
		std::string debug_particle_prefix_;
    std::string log_file_;

    // Maximum phi to compute.
    T phi_max_;

		// These hold particle and grid data for the simulation.
		VectorGridT *grid_velocity_; 
		VectorGridT *grid_velocity_update_;
		VectorGridT *grid_weight_;
		ScalarGridT *grid_phi_;
		ScalarGridT *grid_divergence_;
		ScalarGridT *grid_pressure_;
		ScalarGridInt *grid_marker_;
		ScalarGridBool *grid_source_;
		ParticlesT *particles_;

    // Fluid distribution, store profiled information.
    common::ProfileSliceLocal *profile_;

    // Sources
    std::vector< std::function<void(void)> > sources_;

		// Current frame number to simulate/save.
		int frame_;
		int debug_step_;

    // Intermediate buffer stuff.
    const int log2dim_  = N2;
    const int leaf_dim_ = 1 << N2;
    const int buf_dim_  = leaf_dim_ + 2;
    const int buf_stride_ = buf_dim_ * buf_dim_ * buf_dim_;
    const int buf_num_    = 3*buf_stride_;
    inline int buf_offset(int i, int j, int k) const {
      return (i * buf_dim_ + j) * buf_dim_ + k;
    }

		// Get barycentric weights to transfer particle quantities to a mac grid.
		void GetCoord(const common::Vec3<T> &position, common::Coord &c) const;
		void GetCoordAndWeights(
        const common::Vec3<T> &position, int dim,
				common::Coord &c, common::Vec2<T> weights[3]) const;
    // Initialize particles using initialized phi/marker.
    void InitializePhiAndParticles(std::function<T(const common::Vec3<T>&)>);

    // Compute phi by sweeping.
    void SweepPhiSingle(const int start[3], const int end[3],
                        const int delta[3]);
    T SolveForPhi(T px, T py, T pz, T pc);

    void SetUpOutputAndLoggingFiles();

  public:
    // Accessors for setting, getting.
    void set_grid_velocity(VectorGridT *grid_velocity) {
      grid_velocity_ = grid_velocity;
    }
    void set_grid_velocity_update(VectorGridT *grid_velocity_update) {
      grid_velocity_update_ = grid_velocity_update;
    }
    void set_grid_weight(VectorGridT *grid_weight) {
      grid_weight_ = grid_weight;
    }
    void set_grid_phi(ScalarGridT *grid_phi) {
      grid_phi_ = grid_phi;
    }
    void set_grid_divergence(ScalarGridT *grid_divergence) {
      grid_divergence_ = grid_divergence;
    }
    void set_grid_pressure(ScalarGridT *grid_pressure) {
      grid_pressure_ = grid_pressure;
    }
    void set_grid_marker(ScalarGridInt *grid_marker) {
      grid_marker_ = grid_marker;
    }
    void set_grid_source(ScalarGridBool *grid_source) {
      grid_source_ = grid_source;
    }
    void set_particles(ParticlesT *particles) {
      particles_ = particles;
    }
    void set_profile(common::ProfileSliceLocal *profile) {
      profile_ = profile;
    }
    VectorGridT* grid_velocity() {
      return grid_velocity_;
    }
    VectorGridT* grid_velocity_update() {
      return grid_velocity_update_;
    }
    VectorGridT* grid_weight() {
      return grid_weight_;
    }
    ScalarGridT* grid_phi() {
      return grid_phi_;
    }
    ScalarGridT* grid_divergence() {
      return grid_divergence_;
    }
    ScalarGridT* grid_pressure() {
      return grid_pressure_;
    }
    ScalarGridInt* grid_marker() {
      return grid_marker_;
    }
    ScalarGridBool* grid_source() {
      return grid_source_;
    }
    ParticlesT* particles() {
      return particles_;
    }
    common::ProfileSliceLocal* profile() {
      return profile_;
    }

    // Particle functions to transfer/update single particle.
    void TransferSingleParticleToGrid(
      const ParticleData &pd, const common::Coord c, int total,
      const common::Coord lc, T *v_buffer, T *w_buffer,
      char *mask_buffer) const;
    // Get weights for transferring to grid.
    template<int D>
    void ScatterWeight(const ParticleData &pd, const common::Coord c,
                       const common::Coord lc, T *v_buffer, T *w_buffer,
                       char *mask_buffer) const;
    void UpdateSingleParticleFromGrid(
      ParticleData &pd, const common::Coord c, int total,
      const common::Coord lc, T *v_buffer, T *dv_buffer);
    common::Vec3<T> Interpolate(
      const common::Vec3<T> &pos, const common::Coord c,
      const common::Coord lc, T *v_buffer) const;
    template<int D>
    T InterpolateComponent(
      const common::Vec3<T> &pos, const common::Coord c,
      const common::Coord lc, T *v_buffer) const;

};  // class FlipSimulation

}  // namespace application

#endif  // PROJECTS_FLIP_COMMON_FLIP_SIMULATION_H
