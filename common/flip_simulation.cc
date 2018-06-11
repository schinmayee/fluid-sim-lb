#include <algorithm>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <string>

#include "common/flip_simulation.h"
#include "common/primitives.h"

#define THRESH 1e-5
#define PHI_THRESH_POS 0.5

namespace application {

template<common::Index N1, common::Index N2, typename T>
FlipSimulation<N1, N2, T>::FlipSimulation() : metagrid_() {
  voxel_len_ = 1;
  num_particles_per_cell_ = 8;
  gravity_ = 9.8;
  flip_factor_ = 1;
  debug_ = false;
  frame_ = 0;
  grid_velocity_ = nullptr;
  grid_velocity_update_ = nullptr;
  grid_weight_ = nullptr;
  grid_phi_ = nullptr;
  grid_divergence_ = nullptr;
  grid_marker_ = nullptr;
  grid_solid_ = nullptr;
  particles_ = nullptr;
  profile_ = nullptr;
  phi_max_ = 5*voxel_len_;
}  // FlipSimulation

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SetUpOutputAndLoggingFiles() {
  if (output_dir_ == "") {
    return;
  }
  // Set files and directories.
  if (!boost::filesystem::is_directory(output_dir_)) {
    if (!boost::filesystem::create_directories(output_dir_)) {
      std::cerr << "Could not create output directory " << output_dir_;
      exit(-1);
    }
  }  // if is_directory
  debug_output_dir_ = (boost::filesystem::path(output_dir_) /
                       boost::filesystem::path("debug")).string();
  if (!boost::filesystem::is_directory(debug_output_dir_)) {
    if (!boost::filesystem::create_directories(debug_output_dir_)) {
      std::cerr << "Could not create output directory " << debug_output_dir_;
      exit(-1);
    }
  }  // if is_directory
  // Save data with these prefixes followed by "%08d" formatting
  // for frame numbers.
  grid_velocity_prefix_ = (boost::filesystem::path(output_dir_) /
                           boost::filesystem::path("grid_velocity_")).string();
  grid_phi_prefix_ = (boost::filesystem::path(output_dir_) /
                      boost::filesystem::path("grid_phi_")).string();
  particle_prefix_ = (boost::filesystem::path(output_dir_) /
                      boost::filesystem::path("particles_")).string();
  debug_grid_velocity_prefix_ =
    (boost::filesystem::path(debug_output_dir_) /
     boost::filesystem::path("grid_velocity_")).string();
  debug_grid_phi_prefix_ =
    (boost::filesystem::path(debug_output_dir_) /
     boost::filesystem::path("grid_phi_")).string();
  debug_particle_prefix_ =
    (boost::filesystem::path(debug_output_dir_) /
     boost::filesystem::path("particles_")).string();
  log_file_ = (boost::filesystem::path(output_dir_) /
               boost::filesystem::path("log.txt")).string();
  if (boost::filesystem::exists(log_file_)) {
    boost::filesystem::remove(log_file_);
  }
}  // SetUpOutputAndLoggingFiles

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializeMetadata(
	common::Coord global_dims, common::Coord partitions,
  int rank, T voxel_len, int num_particles_per_cell,
  T gravity, T flip_factor, std::string output_dir, bool debug_flag) {
  metagrid_ = common::MetaGrid(global_dims, partitions, rank);
  VLOG(1) << "Metagrid information:";
  VLOG(1) << metagrid_.string();
  voxel_len_ = voxel_len;
  num_particles_per_cell_ = num_particles_per_cell;
  gravity_ = gravity;
  flip_factor_ = flip_factor;
  output_dir_ = output_dir + "_" + std::to_string(rank);
  debug_ = debug_flag;
  frame_ = 0;
  debug_step_ = 0;
  assert(abs(voxel_len_- 1.0) < 1e-2);
  SetUpOutputAndLoggingFiles();
}  // InitializeMetadata

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializeVariables() {
  CHECK(grid_velocity_ != nullptr);
  CHECK(grid_velocity_update_ != nullptr);
  CHECK(grid_weight_ != nullptr);
  CHECK(grid_phi_ != nullptr);
  CHECK(grid_divergence_ != nullptr);
  CHECK(grid_marker_ != nullptr);
  CHECK(grid_solid_ != nullptr);
  CHECK(particles_ != nullptr);
  grid_velocity_->InitializePartition(metagrid_);
  grid_velocity_update_->InitializePartition(metagrid_);
  grid_weight_->InitializePartition(metagrid_);
  grid_phi_->InitializePartition(metagrid_, 1, phi_max_);
  grid_divergence_->InitializePartition(metagrid_);
  grid_pressure_->InitializePartition(metagrid_);
  grid_marker_->InitializePartition(metagrid_);
  grid_solid_->InitializePartition(metagrid_);
  particles_->InitializePartition(metagrid_);
}  // InitializeVariables

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializeProfiling(common::Coord main_partitions) {
  assert(profile_ != nullptr);
  common::Coord coarse_partitions = metagrid_.partitions();
  common::Coord pidx = metagrid_.partition_id();
  common::Coord start, end;
  for (int i = 0; i < 3; ++i) {
    int width = main_partitions[i]/coarse_partitions[i];
    start[i] = width * pidx[i];
    end[i] = start[i] + width - 1;
  }  // for i
  profile_->Initialize(common::CoordBBox(start, end));
}  // InitializeProfiling

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializePhiAndParticles(
  std::function<T(const common::Vec3<T>&)> ComputePhi,
  const common::Vec3<T> velocity) {
  common::Coord start = metagrid_.start();
  common::Coord end = metagrid_.end();
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        if (grid_solid_->get(i,j,k) == 1) {
          grid_marker_->set(i,j,k, common::SOLID);
          continue;
        }
        common::Vec3<T> cell_origin(
          voxel_len_*T(i), voxel_len_*T(j), voxel_len_*T(k));
        common::Vec3<T> cell_center = cell_origin + T(0.5*voxel_len_);
        T phi = ComputePhi(cell_center);
        if (phi < PHI_THRESH_POS) {
          grid_phi_->set(i,j,k, phi);
          if (phi < 0) {
            grid_marker_->set(i,j,k, common::FLUID);
            if (particles_->num_particles(i, j, k) >= num_particles_per_cell_) {
              continue;
            }
            const int max_attempts = 32;
            int pi = 0;
            int num_attempts = 0;
            for (int num_attempts = 0;
                  num_attempts < max_attempts && pi < num_particles_per_cell_;
                 ++num_attempts) {
              T px = cell_origin.x() + (T(random())/T(RAND_MAX)) * voxel_len_;
              T py = cell_origin.y() + (T(random())/T(RAND_MAX)) * voxel_len_;
              T pz = cell_origin.z() + (T(random())/T(RAND_MAX)) * voxel_len_;
              common::Vec3<T> ppos(px, py, pz);
              if (ComputePhi(ppos) > 0) {
                continue;
              } else {
                ParticleData pd =
                  { .position = ppos, .velocity = velocity,
                    .grid_velocity = velocity };
                particles_->AddParticleToData(i, j, k, pd);
                ++pi;
                ++num_attempts;
              }  // if ... > 0
            }  // for num_attempts
          } // if phi < 0
        }  // if phi < PHI_THRESH_POS
      }  // for k
    }  // for i
  }  // for j
}  // InitializePhiAndParticles

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializeReservoir(T height) {
  std::function<T(const common::Vec3<T>&)> ComputePhi =
    [height](const common::Vec3<T> &position) {
       return(position.y() - height);
     };  // ComputePhi
  InitializePhiAndParticles(ComputePhi, common::Vec3<T>(0));
}  // InitializeReservoir

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializeMultipleWaterDrops(
  T height,
  const std::vector<common::Vec3<T>> &centers,
  const std::vector<T> &radii) {
  VLOG(1) << "Multiple drop initialization for";
  VLOG(1) << metagrid_.string();
  CHECK(centers.size() == radii.size());
  std::function<T(const common::Vec3<T>&)> ComputePhi =
    [height, centers, radii](const common::Vec3<T> &position) {
    T phi = position.y() - height;
    int num_spheres = centers.size();
    for (int i = 0; i < num_spheres; ++i) {
      common::Vec3<T> delta = position - centers[i];
      T phi_sphere = delta.length() - radii[i];
      phi = fmin(phi, phi_sphere);
    }
    return phi;
  };
  InitializePhiAndParticles(ComputePhi, common::Vec3<T>(0));
}  // InitializeMultipleDrops

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializeWaterCuboid(
    const common::Cube<T> &box, const common::Vec3<T> velocity) {
  std::function<T(const common::Vec3<T>&)> ComputePhi =
    [box](const common::Vec3<T> &position) {
      if (box.IsInside(position)) {
        return -1;
      } else {
        return 1;
      }
    };
  InitializePhiAndParticles(ComputePhi, velocity);
}  // InitializeWaterCuboid

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::InitializeWaterSphere(
    const common::Vec3<T> &center, const T radius,
    const common::Vec3<T> velocity) {
  std::function<T(const common::Vec3<T>&)> ComputePhi =
      [center, radius](const common::Vec3<T> &position) {
    common::Vec3<T> delta = position - center;
    T phi = delta.length() - radius;
    return phi;
  };
  InitializePhiAndParticles(ComputePhi, velocity);
}  // InitializeWaterSphere

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SetForces(const common::Vec3<T> &f) {
  fx_ = f.x();
  fy_ = f.y();
  fz_ = f.z();
}  // SetForces

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SetBoxBoundary() {
  common::Coord start = metagrid_.start();
  common::Coord end = metagrid_.end();
  common::Coord global_dims = metagrid_.global_dims();
  // x component, lower
  if (start.x() == 0) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
    for (int k = start[2]-1; k <= end[2]+1; ++k) {
      grid_solid_->set(-1,j,k, 1);
    }  // for k
    }  // for j
  }  // if start.x() == 0
  // x component, upper
  if (end.x()+1 == global_dims[0]) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
    for (int k = start[2]-1; k <= end[2]+1; ++k) {
      grid_solid_->set(global_dims[0],j,k, 1);
    }  // for k
    }  // for j
  }  // if end.x()+1 == global_dims[0]
  // y component, lower
  if (start.y() == 0) {
    for (int i = start[0]-1; i <= end[0]+1; ++i) {
    for (int k = start[2]-1; k <= end[2]+1; ++k) {
      grid_solid_->set(i,-1,k, 1);
    }  // for k
    }  // for i
  }  // if start.y() == 0
  // y component, upper
  if (end.y()+1 == global_dims[1]) {
    for (int i = start[0]-1; i <= end[0]+1; ++i) {
    for (int k = start[2]-1; k <= end[2]+1; ++k) {
      grid_solid_->set(i,global_dims[1],k, 1);
    }  // for k
    }  // for i
  }  // if end.y()+1 == global_dims[1]
  // z component, lower
  if (start.z() == 0) {
    for (int i = start[0]-1; i <= end[0]+1; ++i) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
      grid_solid_->set(i,j,-1, 1);
    }  // for k
    }  // for i
  }  // if start.z() == 0
  // z component, upper
  if (end.z()+1 == global_dims[2]) {
    for (int i = start[0]-1; i <= end[0]+1; ++i) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
      grid_solid_->set(i,j,global_dims[2], 1);
    }  // for k
    }  // for i
  }  // if end.z()+1 == 0
}  // SetBoxBoundary

// Only stationary solids supported right now --- grid_solid_ stores the mask.
// This mask does not change with time.
// This mask is used to enforce that normal velocity at solid interfaces is
// zero, and used to set grid_marker_ which is used for enforcing
// incompressiblity.
template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::AddSolidCuboid(const common::CoordBBox &box) {
  common::CoordBBox domain = metagrid_.bbox();
  // Expand the domain by 1 to capture solid neighbors for fluid in this
  // region.
  domain.expand(1);
  // Get intersection with solid cuboid and set that region.
  domain.intersect(box); 
  grid_solid_->fill(domain, 1, true);
}  // AddSolidCuboid

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::AddWaterSource(const Source<T> &source) {
  sources_.push_back(source);
}  // AddWaterSource

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SetTankBoundary() {
  common::Coord start = metagrid_.start();
  common::Coord end = metagrid_.end();
  common::Coord global_dims = metagrid_.global_dims();
  // x component, lower
  if (start.x() == 0) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
    for (int k = start[2]-1; k <= end[2]+1; ++k) {
      grid_solid_->set(-1,j,k, 1);
    }  // for k
    }  // for j
  }  // if start.x() == 0
  // x component, upper
  if (end.x()+1 == global_dims[0]) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
    for (int k = start[2]-1; k <= end[2]+1; ++k) {
      grid_solid_->set(global_dims[0],j,k, 1);
    }  // for k
    }  // for j
  }  // if end.x()+1 == global_dims[0]
  // y component, lower
  if (start.y() == 0) {
    for (int i = start[0]-1; i <= end[0]+1; ++i) {
    for (int k = start[2]-1; k <= end[2]+1; ++k) {
      grid_solid_->set(i,-1,k, 1);
    }  // for k
    }  // for i
  }  // if start.y() == 0
  // No solid for upper y as it is an open tank
  // z component, lower
  if (start.z() == 0) {
    for (int i = start[0]-1; i <= end[0]+1; ++i) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
      grid_solid_->set(i,j,-1, 1);
    }  // for k
    }  // for i
  }  // if start.z() == 0
  // z component, upper
  if (end.z()+1 == global_dims[2]) {
    for (int i = start[0]-1; i <= end[0]+1; ++i) {
    for (int j = start[1]-1; j <= end[1]+1; ++j) {
      grid_solid_->set(i,j,global_dims[2], 1);
    }  // for k
    }  // for i
  }  // if end.z()+1 == 0
}  // SetTankBoundary

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ClearGridData() {
  // Do not clear solids and sources.
  grid_velocity_->Clear(); 
  grid_velocity_update_->Clear();
  grid_weight_->Clear();
  grid_phi_->Clear();
  grid_divergence_->Clear();
  grid_marker_->Clear();
  grid_pressure_->Clear();
}  // ClearGridData

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::Profile(
  float compute_time, float scatter_gather_time, bool use_geometric) {
  CHECK(profile_ != nullptr);
  CHECK(grid_marker_ != nullptr);
  profile_->Clear();
  // Scaling factor.
  const common::Coord global_dims = metagrid_.global_dims();
  float scale = 1.0;
  for (int i = 0; i < 3; ++i) {
    scale = scale * (T(global_dims[i])/128.0);
  }
  scale = 1.0/scale;
  LOG(INFO) << "Global dims " << common::ToString(global_dims) <<
               " and scaling " << scale;
  // Compute partitions that this region contains.
  common::CoordBBox partitions = profile_->get_local_partitions();
  common::Coord pwidth = partitions.dim();
  common::Coord pstart = partitions.min();
  common::Coord pend   = partitions.max();
  // Compute dimensions of each partition to profile.
  common::Coord dims = metagrid_.local_dims();
  common::Coord width(0);
  for (int d = 0; d < 3; ++d) {
    assert(pwidth[d] > 0);
    width[d] = dims[d]/pwidth[d];
  }
  // Allocate a new slice.
  profile_->set_compute_time(compute_time);
  profile_->set_scatter_gather_time(scatter_gather_time);
  // Now count and store number of fluid cells.
  for (int i = pstart[0]; i <= pend[0]; ++i) {
    for (int j = pstart[1]; j <= pend[1]; ++j) {
      for (int k = pstart[2]; k <= pend[2]; ++k) {
        if (use_geometric) {
          profile_->set_count(i, j, k, 0);
        } else {
          common::Coord start(i*width[0], j*width[1], k*width[2]);
          common::Coord end = start + width - common::Coord(1,1,1);
          common::CoordBBox box(start, end);
          int num_fluid_cells = grid_marker_->Count(box, common::FLUID);
          int scaled_cells = scale * T(num_fluid_cells);
          profile_->set_count(i, j, k, scaled_cells);
        }
      }  // for k
    }  // for j
  }  // for i
}  // Profile

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::PrintProfilingData(
  std::string prefix, bool print_rank) {
  if (profile_ == nullptr) {
    return;
  }
  LOG(INFO) << "Profiling data for " << prefix;
  if (print_rank) {
    LOG(INFO) << profile_->Print(&metagrid_);
  } else {
    LOG(INFO) << profile_->Print(nullptr);
  }
}  // PrintProfilingData

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ClearProfiledData() {
  profile_->Clear();
}  // ClearProfiledData

template<common::Index N1, common::Index N2, typename T>
T FlipSimulation<N1, N2, T>::CFL() {
  common::Coord start = metagrid_.start();
  common::Coord end = metagrid_.end();
  T v_inf = 0.0;
  for (typename VectorGridT::GridT::ValueOnCIter iter =
	  grid_velocity_->data()->cbeginValueOn(); iter.test(); ++iter) {
    common::Coord c = iter.getCoord();
    // Make sure this point is in bounds for each component, and then take
    // the maximum absolute -- ignore boundary/ghost values and staggered grid
    // values outside range.
    if (c[0] < start[0] || c[0] > end[0] ||
        c[1] < start[1] || c[1] > end[1] ||
        c[2] < start[2] || c[2] > end[2]) {
      continue;
    }
    if (c[1] != end[1] && c[2] != end[2]) {
      v_inf = fmax(v_inf, fabs((iter.getValue())[0]));
    }
    if (c[0] != end[0] && c[2] != end[2]) {
      v_inf = fmax(v_inf, fabs((iter.getValue())[1]));
    }
    if (c[0]!= end[0] && c[1] != end[1]) {
      v_inf = fmax(v_inf, fabs((iter.getValue())[2]));
    }
  }  // for iter
  T v = fmax(sqrt(0.5*gravity_), v_inf);
  if (v < 1e-2) v = 1e-2;
  return 1.0/v;
}  // CFL

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SaveGridVelocities() {
  for (typename VectorGridT::GridT::ValueOnCIter iter =
	  grid_velocity_->data()->cbeginValueOn(); iter.test(); ++iter) {
	grid_velocity_update_->set(iter.getCoord(), iter.getValue());
  }  // for iter
}  // SaveGridVelocities

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::MoveParticles(T dt, bool defragment) {
  particles_->template StepInGrid<VectorGridT, true, 1, 1>(
    *grid_velocity_, dt);
  particles_->EvaluateActiveBoundingBox();
}  // MoveParticles

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::AddGravity(T dt) {
  for (typename VectorGridT::GridT::ValueOnIter iter =
    grid_velocity_->data()->beginValueOn(); iter.test(); ++iter) {
    common::Vec3<T> v = iter.getValue();
    v[0] += dt * fx_;
    v[1] -= dt * fy_;
    v[2] += dt * fz_;
    iter.setValue(v);
  }
}  // AddGravity

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::AdvectPhi() {
  // TODO: later, at the end, after basic simulation works
}

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ReinitializePhi() {
  std::cerr << "ReinitializePhi unimplemented!";
  exit(-1);
}  // ReinitializePhi

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ReinitializePhiSimple() {
  // TODO: later
  // First set phi in fluid to -voxel_len/2.
  T phi_fluid = -0.5;
  for (typename ScalarGridInt::GridT::ValueOnCIter iter =
      grid_marker_->data()->cbeginValueOn(); iter.test(); ++iter) {
    if (iter.getValue() == common::FLUID) {
      grid_phi_->set(iter.getCoord(), phi_fluid);
    }
  }  // for iter
}  // ReinitializePhiOutside

// Solution for this should be (after solving quadratic representing sum of
// gradient):
// s = sqrt(3*voxel_len_^2 - 2(px^2 + py^2 + pz^2 - px*py - py*pz - px*pz))
// phi = ((px + py + pz) +- s)/3.
template<common::Index N1, common::Index N2, typename T>
T FlipSimulation<N1, N2, T>::SolveForPhi(T px, T py, T pz, T pc) {
  T phi = fmin(fmin(px, py), pz) + 1;
  if (phi > fmax(fmax(px, py), pz)) {
    T s = 3 - 2*(px*(px-py) + py*(py-pz) + pz*(pz-px));
    assert(s >= 0);
    s = sqrt(s);
    phi = ((px + py + pz) + s)/3;
    assert(phi >= px || phi >= py || phi >= pz);
  }
  return fmin(phi, pc);
}  // SolveForPhi

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SweepPhiSingle(
  const int start[3], const int end[3], const int delta[3]) {
  assert(abs(start[0]-end[0])+1 == metagrid_.dims(0));
  assert(abs(start[1]-end[1])+1 == metagrid_.dims(1));
  assert(abs(start[2]-end[2])+1 == metagrid_.dims(2));
  for (int i = start[0]; i != end[0]; i += delta[0]) {
    for (int j = start[1]; j != end[1]; j += delta[1]) {
      for (int k = start[2]; k != end[2]; k += delta[2]) {
        if (grid_marker_->get(i,j,k) != common::FLUID) {
          T phi_x = grid_phi_->get(i-delta[0],j,k);
          T phi_y = grid_phi_->get(i,j-delta[1],k);
          T phi_z = grid_phi_->get(i,j,k-delta[2]);
          if (phi_x >= phi_max_ && phi_y >= phi_max_ && phi_z >= phi_max_) {
            // There is no need to compute phi for this case since it is going
            // to be greater than phi_max_.
            continue;
          }
          T phi = SolveForPhi(phi_x, phi_y, phi_z, grid_phi_->get(i,j,k));
          if (phi < phi_max_) {
            grid_phi_->set(i,j,k,phi);
          }  // phi < phi_upper_thresh
        }  // if common::FLUID
      }  // for k
    }  // for j
  }  // for i
}  // SweepPhiSingle

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SweepPhi() {
  // Compute phi by sweeping in all 8 directions.
  common::Coord lo = metagrid_.start();
  common::Coord hi = metagrid_.end();
  {
    int start[3] = {lo[0], lo[1], lo[2]};
    int end[3] =   {hi[0], hi[1], hi[2]};
    int delta[3] = {1,1,1};
    SweepPhiSingle(start, end, delta);
  }
  {
    int start[3] = {lo[0], lo[1], hi[2]};
    int end[3] =   {hi[0], hi[1], lo[2]};
    int delta[3] = {1,1,-1};
    SweepPhiSingle(start, end, delta);
  }
  {
    int start[3] = {lo[0], hi[1], lo[2]};
    int end[3] =   {hi[0], lo[1], hi[2]};
    int delta[3] = {1,-1,1};
    SweepPhiSingle(start, end, delta);
  }
  {
    int start[3] = {lo[0], hi[1], hi[2]};
    int end[3] =   {hi[0], lo[1], lo[2]};
    int delta[3] = {1,-1,-1};
    SweepPhiSingle(start, end, delta);
  }
  {
    int start[3] = {hi[0], lo[1], lo[2]};
    int end[3] =   {lo[0], hi[1], hi[2]};
    int delta[3] = {-1,1,1};
    SweepPhiSingle(start, end, delta);
  }
  {
    int start[3] = {hi[0], lo[1], hi[2]};
    int end[3] =   {lo[0], hi[1], lo[2]};
    int delta[3] = {-1,1,-1};
    SweepPhiSingle(start, end, delta);
  }
  {
    int start[3] = {hi[0], hi[1], lo[2]};
    int end[3] =   {lo[0], lo[1], hi[2]};
    int delta[3] = {-1,-1,1};
    SweepPhiSingle(start, end, delta);
  }
  {
    int start[3] = {hi[0], hi[1], hi[2]};
    int end[3] =   {lo[0], lo[1], lo[2]};
    int delta[3] = {-1,-1,-1};
    SweepPhiSingle(start, end, delta);
  }
}  // SweepPhi

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ExtendVelocity() {
  // TODO: later, after ReinitializePhi
}

#define SET_COMPONENT_TO_ZERO(coord, dim) (  \
{  \
  if (grid_velocity_->isOn(coord)) {\
    common::Vec3<T> v = grid_velocity_->get(coord);  \
    v[dim] = 0;  \
    grid_velocity_->set(coord, v);  \
  }  \
}  \
)

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ApplyBoundaryConditions() {
  // Set solid marker, and velocity where there is solid to zero.
  for (typename ScalarGridInt::GridT::ValueOnCIter iter =
       grid_solid_->data()->cbeginValueOn(); iter.test(); ++iter) {
    CHECK_EQ(iter.getValue(), 1);
    const common::Coord c = iter.getCoord();
    grid_marker_->set(c, common::SOLID);
    // x lower
    if (grid_solid_->get(c.x()-1, c.y(), c.z()) == 0) {
      if (grid_velocity_->isOn(c)) {
        common::Vec3<T> v = grid_velocity_->get(c);
        v[0] = 0;
        grid_velocity_->set(c, v);
      }
    }
    // y lower
    if (grid_solid_->get(c.x(), c.y()-1, c.z()) == 0) {
      if (grid_velocity_->isOn(c)) {
        common::Vec3<T> v = grid_velocity_->get(c);
        v[1] = 0;
        grid_velocity_->set(c, v);
      }
    }
    // z lower
    if (grid_solid_->get(c.x(), c.y(), c.z()-1) == 0) {
      if (grid_velocity_->isOn(c)) {
        common::Vec3<T> v = grid_velocity_->get(c);
        v[2] = 0;
        grid_velocity_->set(c, v);
      }
    }
    // x upper
    if (grid_solid_->get(c.x()+1, c.y(), c.z()) == 0) {
      if (grid_velocity_->isOn(c.x()+1, c.y(), c.z())) {
        common::Vec3<T> v = grid_velocity_->get(c.x()+1, c.y(), c.z());
        v[0] = 0;
        grid_velocity_->set(c.x()+1, c.y(), c.z(), v);
      }
    }
    // y upper
    if (grid_solid_->get(c.x(), c.y()+1, c.z()) == 0) {
      if (grid_velocity_->isOn(c.x(), c.y()+1, c.z())) {
        common::Vec3<T> v = grid_velocity_->get(c.x(), c.y()+1, c.z());
        v[1] = 0;
        grid_velocity_->set(c.x(), c.y()+1, c.z(), v);
      }
    }
    // z upper
    if (grid_solid_->get(c.x(), c.y(), c.z()+1) == 0) {
      if (grid_velocity_->isOn(c.x(), c.y(), c.z()+1)) {
        common::Vec3<T> v = grid_velocity_->get(c.x(), c.y(), c.z()+1);
        v[2] = 0;
        grid_velocity_->set(c.x(), c.y(), c.z()+1, v);
      }
    }
  }  // for iter ... grid_solid_
}  // ApplyBoundaryConditions

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ComputeDivergence() {
  common::Coord start = metagrid_.start();
  common::Coord end = metagrid_.end();
  // First compute divergence (rhs). All differentials are in index space.
  for (int i = start[0]; i <= end[0]; ++i) {
		for (int j = start[1]; j <= end[1]; ++j) {
			for (int k = start[2]; k <= end[2]; ++k) {
				if (grid_marker_->get(i,j,k) == common::FLUID) {
					T divergence = 0;
					divergence += grid_velocity_->get(i+1,j,k)[0] -
                        grid_velocity_->get(i,j,k)[0];
					divergence += grid_velocity_->get(i,j+1,k)[1] -
                        grid_velocity_->get(i,j,k)[1];
					divergence += grid_velocity_->get(i,j,k+1)[2] -
                        grid_velocity_->get(i,j,k)[2];
					grid_divergence_->set(i, j, k, divergence);
				}  // if common::FLUID
			}  // for k
		}  // for j
  }  // for i
}  // ComputeDivergence

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::UpdateVelocityFromPressure() {
  // Update velocity by taking gradient of computed pressure.
  common::Coord start = metagrid_.start();
  common::Coord end = metagrid_.end();
	for (int i = start[0]; i <= end[0]+1; ++i) {
		for (int j = start[1]; j <= end[1]+1; ++j) {
			for (int k = start[2]; k <= end[2]+1; ++k) {
				common::Vec3<T> v = grid_velocity_->get(i,j,k);
        if (grid_marker_->get(i,j,k) == common::SOLID) {
          continue;
        }
        bool update = false;
				if ((grid_marker_->get(i-1,j,k) == common::FLUID ||
						 grid_marker_->get(i,j,k) == common::FLUID) &&
            grid_marker_->get(i-1,j,k) != common::SOLID) {
					v[0] += (grid_pressure_->get(common::Coord(i,j,k)) -
                   grid_pressure_->get(common::Coord(i-1,j,k)));
          update = true;
				}
				if ((grid_marker_->get(i,j-1,k) == common::FLUID ||
						 grid_marker_->get(i,j,k) == common::FLUID) &&
            grid_marker_->get(i,j-1,k) != common::SOLID) {
					v[1] += (grid_pressure_->get(common::Coord(i,j,k)) -
                   grid_pressure_->get(common::Coord(i,j-1,k)));
          update = true;
				}
				if ((grid_marker_->get(i,j,k-1) == common::FLUID ||
						 grid_marker_->get(i,j,k) == common::FLUID) &&
            grid_marker_->get(i,j,k-1) != common::SOLID) {
					v[2] += (grid_pressure_->get(common::Coord(i,j,k)) -
                   grid_pressure_->get(common::Coord(i,j,k-1)));
          update = true;
				}
        if (update) {
          grid_velocity_->set(i,j,k,v);
        }
			}  // for k
		}  // for j
	}  // for i
}  // UpdateVelocityFromPressure

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::GetVelocityUpdate() {
  const common::Coord start = metagrid_.start();
  const common::Coord end = metagrid_.end();
	for (int i = start[0]; i <= end[0]+1; ++i) {
		for (int j = start[1]; j <= end[1]+1; ++j) {
			for (int k = start[2]; k <= end[2]+1; ++k) {
        const common::Coord ci(i,j,k);
        if (grid_velocity_->isOn(ci) && !grid_velocity_update_->isOn(ci)) {
          CHECK(false) << common::ToString(ci) << " "
                       << common::ToString(start) << " "
                       << common::ToString(end);
        }
        if (grid_velocity_->isOn(ci) && grid_velocity_update_->isOn(ci)) {
          const common::Vec3<T> update =
            grid_velocity_->get(ci) - grid_velocity_update_->get(ci);
          grid_velocity_update_->set(ci, update);
        }
      }  // for k
    }  // for j
  }  // for i
}  // GetVelocityUpdate

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ExtendGridVelocityIntoBoundaryRegion() {
  // TODO: later, to prevent stuck particles at boundary region
}  // ExtendGridVelocityIntoBoundaryRegion

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::TransferParticlesToGrid() {
  // Get weights and accumulate velocity for each component (dimension).
  // Transfer weighted velocity from each particle to grid.
  //
  common::Coord start = metagrid_.start();
  common::Coord end = metagrid_.end();

  typedef typename ParticlesT::Address Address;

  typename ParticlesT::GridT::Ptr data = particles_->data();
  int total = 0;

  // d, then x,y,z
  T *v_buffer = new T[buf_num_];
  T *w_buffer = new T[buf_num_];
  char *mask_buffer = new char[buf_num_];

  for (auto leaf_iter = data->tree().cbeginLeaf(); leaf_iter.test(); ++leaf_iter) {

    const common::Coord origin = leaf_iter->origin() - common::Coord(1,1,1);

    // Reset buffer for every leaf as buffer represents a leaf.
    for (int i = 0; i < buf_num_; ++i) {
      v_buffer[i] = 0;
      w_buffer[i] = 0;
      mask_buffer[i] = 0;
    }

    for (auto iter = leaf_iter->cbeginValueOn(); iter.test(); ++iter) {
      Address addr = iter.getValue();
      if (addr == 0) {
        continue;
      }  // if addr == 0

      typename ParticlesT::ParticleList *list =
        reinterpret_cast<typename ParticlesT::ParticleList*>(addr);
      if (list->size() == 0) {
        continue;
      }

      const common::Coord c = iter.getCoord();
      const common::Coord lc = c - origin;

      // Mark fluid.
      if (c[0] >= start[0] && c[0] <= end[0] &&
          c[1] >= start[1] && c[1] <= end[1] &&
          c[2] >= start[2] && c[2] <= end[2]) {
        grid_marker_->set(c, common::FLUID);
      }

      // Now scatter weight and weighted velocity into buffers.
      std::function<void(const ParticleData&, const common::Coord, int)>
        transfer_op = std::bind(
          &FlipSimulation<N1, N2, T>::TransferSingleParticleToGrid, this,
          std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
          lc, v_buffer, w_buffer, mask_buffer);
      total = list->ForEachConst(transfer_op, iter.getCoord(), total);

    }  // for iter

    // Save buffer data for the leaf.
    for (int i = 0; i < buf_dim_; ++i) {
      for (int j = 0; j < buf_dim_; ++j) {
        for (int k = 0; k < buf_dim_; ++k) {
          common::Coord to_set = origin + common::Coord(i,j,k);
          int offset = buf_offset(i, j, k);
          const T *v_base = v_buffer + offset;
          const T *w_base = w_buffer + offset;
          if (mask_buffer[offset] == 0 &&
              mask_buffer[offset+buf_stride_] == 0 &&
              mask_buffer[offset+2*buf_stride_] == 0) {
            continue;
          }
          common::Vec3<T> v = grid_velocity_->get(to_set);
          v += common::Vec3<T>(v_base[0], v_base[buf_stride_], v_base[2*buf_stride_]);
          grid_velocity_->set(to_set, v);
          common::Vec3<T> w = grid_weight_->get(to_set);
          w += common::Vec3<T>(w_base[0], w_base[buf_stride_], w_base[2*buf_stride_]);
          grid_weight_->set(to_set, w);
        }  // for k
      }  // for j
    }  // for i

  }  // for leaf_iter

  // Free buffers.
  delete v_buffer;
  delete w_buffer;

  assert(total == particles_->size());

  // Now normalize grid velocity with the weights.
  #pragma GCC_ivdep
  for (typename VectorGridT::GridT::ValueOnIter iter =
	  grid_velocity_->data()->beginValueOn(); iter.test(); ++iter) {
    common::Vec3<T> v = iter.getValue();
    common::Vec3<T> wsum = grid_weight_->get(iter.getCoord());
    if (wsum[0] < THRESH) {
      v[0] = 0;
    } else {
      v[0] /= wsum[0];
    }
    if (wsum[1] < THRESH) {
      v[1] = 0;
    } else {
      v[1] /= wsum[1];
    }
    if (wsum[2] < THRESH) {
      v[2] = 0;
    } else {
      v[2] /= wsum[2];
    }
    iter.setValue(v);
  }  // for iter
}  // TransferParticlesToGrid

template<common::Index N1, common::Index N2, typename T>
template<int D>
void FlipSimulation<N1, N2, T>::ScatterWeight(
  const ParticleData &pd, const common::Coord c,
  const common::Coord lc, T *v_buffer, T *w_buffer,
  char *mask_buffer) const {

  assert(int(floor(pd.position[0])) == c[0]);
  assert(int(floor(pd.position[1])) == c[1]);
  assert(int(floor(pd.position[2])) == c[2]);

  common::Vec3<T> lower_vel_position(c[0], c[1], c[2]);
  int start[3] = {lc[0], lc[1], lc[2]};

  // Determine position where lower velocity resides.
  // Branch is decided at compile time.
  if (D == 0) {
    lower_vel_position += common::Vec3<T>(0, 0.5, 0.5);
    if (pd.position[1] < lower_vel_position[1]) {
      lower_vel_position[1] -= T(1);
      start[1] -= 1;
    }
    if (pd.position[2] < lower_vel_position[2]) {
      lower_vel_position[2] -= T(1);
      start[2] -= 1;
    }
  }
  if (D == 1) {
    lower_vel_position += common::Vec3<T>(0.5, 0, 0.5);
    if (pd.position[0] < lower_vel_position[0]) {
      lower_vel_position[0] -= T(1);
      start[0] -= 1;
    }
    if (pd.position[2] < lower_vel_position[2]) {
      lower_vel_position[2] -= T(1);
      start[2] -= 1;
    }
  }
  if (D == 2) {
    lower_vel_position += common::Vec3<T>(0.5, 0.5, 0);
    if (pd.position[0] < lower_vel_position[0]) {
      lower_vel_position[0] -= T(1);
      start[0] -= 1;
    }
    if (pd.position[1] < lower_vel_position[1]) {
      lower_vel_position[1] -= T(1);
      start[1] -= 1;
    }
  }

  // Compute all weights.
  T weights[2][2][2];
  {
    T x = pd.position.x() - lower_vel_position.x();
    T y = pd.position.y() - lower_vel_position.y();
    T z = pd.position.z() - lower_vel_position.z();
    T xy = x*y, yz = y*z, xz = z*x, xyz = xy*z;
    weights[1][1][1] = xyz;
    weights[1][1][0] = xy - xyz;
    weights[1][0][1] = xz - xyz;
    weights[1][0][0] = x - xz - weights[1][1][0];
    weights[0][1][1] = yz - xyz;
    weights[0][1][0] = y - yz - weights[1][1][0];
    weights[0][0][1] = z - yz - weights[1][0][1];
    weights[0][0][0] = 1 - x - y + xy - weights[0][0][1];
  }

  // Scatter.
  #pragma GCC_ivdep
  for (int di = 0; di < 2; ++di) {
    int i = di + start[0];
    #pragma GCC_ivdep
    for (int dj = 0; dj < 2; ++dj) {
      int j = dj + start[1];
      #pragma GCC_ivdep
      for (int dk = 0; dk < 2; ++dk) {
        int k = dk + start[2];
        int idx = buf_offset(i, j, k);
        v_buffer[idx] += weights[di][dj][dk] * pd.velocity[D];
        w_buffer[idx] += weights[di][dj][dk];
        mask_buffer[idx] = 1;
      }  // for k
    }  // for j
  }  // for i
}  // ScatterWeight

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::TransferSingleParticleToGrid(
  const ParticleData &pd, const common::Coord c, int total,
  const common::Coord lc, T *v_buffer, T *w_buffer,
  char *mask_buffer) const {
  // Scatter each component.
  ScatterWeight<0>(pd, c, lc, v_buffer, w_buffer, mask_buffer);
  ScatterWeight<1>(pd, c, lc, v_buffer+1*buf_stride_, w_buffer+1*buf_stride_,
                   mask_buffer+1*buf_stride_);
  ScatterWeight<2>(pd, c, lc, v_buffer+2*buf_stride_, w_buffer+2*buf_stride_,
                   mask_buffer+2*buf_stride_);
}  // TransferSingleParticleToGrid

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::UpdateParticlesFromGrid() {
  // TODO: add fill boundary cells for interpolating
  // to particles at the boundary.
  // Interpolate update to particles.
  grid_velocity_->ExtrapolateToBoundary();
  grid_velocity_update_->ExtrapolateToBoundary();

  typedef typename ParticlesT::Address Address;

  typename ParticlesT::GridT::Ptr data = particles_->data();
  int total = 0;

  // d, then x,y,z
  T *v_buffer  = new T[buf_num_];
  T *dv_buffer = new T[buf_num_];

  for (auto leaf_iter = data->tree().cbeginLeaf(); leaf_iter.test(); ++leaf_iter) {

    const common::Coord origin = leaf_iter->origin() - common::Coord(1,1,1);

    // Copy buffer data..
    for (int i = 0; i < buf_dim_; ++i) {
      for (int j = 0; j < buf_dim_; ++j) {
        for (int k = 0; k < buf_dim_; ++k) {
          common::Coord to_get = origin + common::Coord(i,j,k);
          int offset = buf_offset(i, j, k);
          T *v_base  = v_buffer + offset;
          T *dv_base = dv_buffer + offset;
          common::Vec3<T> v  = grid_velocity_->get(to_get);
          v_base[0]              = v[0];
          v_base[buf_stride_]    = v[1];
          v_base[2*buf_stride_]  = v[2];
          common::Vec3<T> dv = grid_velocity_update_->get(to_get);
          dv_base[0]             = dv[0];
          dv_base[buf_stride_]   = dv[1];
          dv_base[2*buf_stride_] = dv[2];
        }  // for k
      }  // for j
    }  // for i

    for (auto iter = leaf_iter->cbeginValueOn(); iter.test(); ++iter) {
      Address addr = iter.getValue();
      if (addr == 0) {
        continue;
      }  // if addr == 0

      typename ParticlesT::ParticleList *list =
        reinterpret_cast<typename ParticlesT::ParticleList*>(addr);
      if (list->size() == 0) {
        continue;
      }

      const common::Coord c = iter.getCoord();
      const common::Coord lc = c - origin;

      // Now scatter weight and weighted velocity into buffers.
      std::function<void(ParticleData&, const common::Coord, int)>
        update_op = std::bind(
          &FlipSimulation<N1, N2, T>::UpdateSingleParticleFromGrid, this,
          std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
          lc, v_buffer, dv_buffer);
      total = list->ForEach(update_op, iter.getCoord(), total);

    }  // for iter

  }  // for leaf_iter
}  // UpdateParticlesFromGrid

template<common::Index N1, common::Index N2, typename T>
template<int D>
T FlipSimulation<N1, N2, T>::InterpolateComponent(
  const common::Vec3<T> &pos, const common::Coord c,
  const common::Coord lc, T *v_buffer) const {
  assert(int(floor(pos[0])) == c[0]);
  assert(int(floor(pos[1])) == c[1]);
  assert(int(floor(pos[2])) == c[2]);

  common::Vec3<T> lower_vel_position(c[0], c[1], c[2]);
  int start[3] = {lc[0], lc[1], lc[2]};

  // Determine position where lower velocity resides.
  // Branch is decided at compile time.
  if (D == 0) {
    lower_vel_position += common::Vec3<T>(0, 0.5, 0.5);
    if (pos[1] < lower_vel_position[1]) {
      lower_vel_position[1] -= T(1);
      start[1] -= 1;
    }
    if (pos[2] < lower_vel_position[2]) {
      lower_vel_position[2] -= T(1);
      start[2] -= 1;
    }
  }
  if (D == 1) {
    lower_vel_position += common::Vec3<T>(0.5, 0, 0.5);
    if (pos[0] < lower_vel_position[0]) {
      lower_vel_position[0] -= T(1);
      start[0] -= 1;
    }
    if (pos[2] < lower_vel_position[2]) {
      lower_vel_position[2] -= T(1);
      start[2] -= 1;
    }
  }
  if (D == 2) {
    lower_vel_position += common::Vec3<T>(0.5, 0.5, 0);
    if (pos[0] < lower_vel_position[0]) {
      lower_vel_position[0] -= T(1);
      start[0] -= 1;
    }
    if (pos[1] < lower_vel_position[1]) {
      lower_vel_position[1] -= T(1);
      start[1] -= 1;
    }
  }

  const int i = start[0];
  const int j = start[1];
  const int k = start[2];
  T w[3] = { pos[0] - lower_vel_position[0],
             pos[1] - lower_vel_position[1],
             pos[2] - lower_vel_position[2] };
  assert(w[0] >= 0 && w[0] <= 1);
  assert(w[1] >= 0 && w[1] <= 1);
  assert(w[2] >= 0 && w[2] <= 1);
  const T v000 = v_buffer[buf_offset(i,j,k)];
  const T v001 = v_buffer[buf_offset(i,j,k+1)];
  const T v010 = v_buffer[buf_offset(i,j+1,k)];
  const T v011 = v_buffer[buf_offset(i,j+1,k+1)];
  const T v100 = v_buffer[buf_offset(i+1,j,k)];
  const T v101 = v_buffer[buf_offset(i+1,j,k+1)];
  const T v110 = v_buffer[buf_offset(i+1,j+1,k)];
  const T v111 = v_buffer[buf_offset(i+1,j+1,k+1)];
  T v00 = v000 + w[2]*(v001-v000);
  T v01 = v010 + w[2]*(v011-v010);
  T v10 = v100 + w[2]*(v101-v100);
  T v11 = v110 + w[2]*(v111-v110);
  T v0 = v00 + w[1]*(v01-v00);
  T v1 = v10 + w[1]*(v11-v10);
  T v  = v0 + w[0]*(v1-v0);
  return v;
}

template<common::Index N1, common::Index N2, typename T>
common::Vec3<T> FlipSimulation<N1, N2, T>::Interpolate(
  const common::Vec3<T> &pos, const common::Coord c,
  const common::Coord lc, T *v_buffer) const {
  T v0 = InterpolateComponent<0>(pos, c, lc, v_buffer);
  T v1 = InterpolateComponent<1>(pos, c, lc, v_buffer+1*buf_stride_);
  T v2 = InterpolateComponent<2>(pos, c, lc, v_buffer+2*buf_stride_);
  return common::Vec3<T>(v0, v1, v2);
}

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::UpdateSingleParticleFromGrid(
  ParticleData &pd, const common::Coord c, int total,
  const common::Coord lc, T *v_buffer, T *dv_buffer) {
  common::Vec3<T> pv = pd.velocity;
  common::Vec3<T> v = Interpolate(pd.position, c, lc, v_buffer);
  common::Vec3<T> dv = Interpolate(pd.position, c, lc, dv_buffer);
  //common::Vec3<T> v = grid_velocity_->Interpolate(pd.position, 1);
  //common::Vec3<T> dv = grid_velocity_update_->Interpolate(pd.position, 1);
  pd.velocity = v + flip_factor_ * (dv + pv - v);
  pd.grid_velocity = v;
}  // UpdateSingleParticleFromGrid

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::ResampleParticles() {
  // TODO: later, after basic simulation
}

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::DeleteParticlesInSolids() {
  particles_->DeleteParticlesInSolids(*grid_solid_);
}  // DeleteParticlesInSolids

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::DeleteOutsideParticles() {
  particles_->DeleteOutsideParticles();
}  // DeleteOutsideParticles

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::AddContributionFromSources(T t) {
  for (const auto source : sources_) {
    if (!source.is_active(t)) {
      continue;
    }
    // Mark regions that are covered by sources as fluid and set the velocity.
    const common::Vec3<T> velocity = source.get_velocity(t);
    const common::Cube<T> box = source.get_box(t);
    const common::Vec3<T> lo = box.get_start();
    const common::Vec3<T> hi = box.get_end();
    const common::Coord start_source(lo[0], lo[1], lo[2]);
    const common::Coord end_source(hi[0], hi[1], hi[2]);
    common::CoordBBox domain = metagrid_.bbox();
    domain.intersect(common::CoordBBox(start_source, end_source));
    const common::Coord start = domain.min();
    const common::Coord end = domain.max();
    for (int i = start[0]; i <= end[0]; ++i) {
      for (int j = start[1]; j <= end[1]; ++j) {
        for (int k = start[2]; k <= end[2]; ++k) {
          common::Vec3<T> cell_origin(
            voxel_len_*T(i), voxel_len_*T(j), voxel_len_*T(k));
          // Mark grid as fluid.
          grid_marker_->set(i,j,k, common::FLUID);
          // Set grid velocity.
          {
            // x lower, y lower, z lower
            grid_velocity_->set(i,j,k, velocity);
            // x upper
            common::Vec3<T> vx = grid_velocity_->get(i+1,j,k);
            vx[0] = velocity[0];
            grid_velocity_->set(i+1,j,k, vx);
            // y upper
            common::Vec3<T> vy = grid_velocity_->get(i,j+1,k);
            vy[1] = velocity[1];
            grid_velocity_->set(i,j+1,k, vy);
            // x upper
            common::Vec3<T> vz = grid_velocity_->get(i,j,k+1);
            vz[2] = velocity[2];
            grid_velocity_->set(i,j,k+1, vz);
          }
          // Add particles.
          const int max_attempts = 32;
          int pi = 0;
           for (int num_attempts = 0;
                num_attempts < max_attempts && pi < num_particles_per_cell_;
                ++num_attempts) {
            T px = cell_origin.x() + (T(random())/T(RAND_MAX)) * voxel_len_;
            T py = cell_origin.y() + (T(random())/T(RAND_MAX)) * voxel_len_;
            T pz = cell_origin.z() + (T(random())/T(RAND_MAX)) * voxel_len_;
            common::Vec3<T> ppos(px, py, pz);
            if (!box.IsInside(ppos)) {
              continue;
            }
            ParticleData pd =
              { .position = ppos, .velocity = velocity,
                .grid_velocity = velocity };
            particles_->AddParticleToData(i, j, k, pd);
            ++pi;
          }  // for num_attempts
        }  // for k
      }  // for j
    }  // for i
  }  // for source
}  // AddContributionFromSources

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::DeleteFluidInSources(T t) {
  std::vector< common::Cube<T> > regions;
  for (auto source : sources_) {
    const common::Cube<T> box = source.get_box(t);
    regions.push_back(box);
  }  // for source
  // Only delete particles. Rest of the grid quantities will be reset anyways.
  particles_->DeleteParticlesFromRegions(regions);
}  // DeleteFluidInSources

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SaveDataToFile() {
  VLOG(1) << "Saving frame " << frame_ << " ...";
  char frame_number_c[10];
  sprintf(frame_number_c, "%08d", frame_);
  std::string frame_number(frame_number_c);
  //grid_velocity_->Write(grid_velocity_prefix_ + frame_number);
  //grid_phi_->Write(grid_phi_prefix_ + frame_number);
  particles_->Write(particle_prefix_ + frame_number, false);
}  // SaveDataToFile

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::SaveDebugDataToFile() {
  assert(debug_);
  {
    VLOG(1) << "Saving binary, debug step " << debug_step_;
    char debug_step_c[10];
    sprintf(debug_step_c, "%08d", debug_step_);
    std::string debug_step_str(debug_step_c);
    grid_velocity_->Write(debug_grid_velocity_prefix_ + debug_step_str);
    grid_phi_->Write(debug_grid_phi_prefix_ + debug_step_str);
    particles_->Write(debug_particle_prefix_ + debug_step_str);
  }
  {
    VLOG(1) << "Saving text, debug step " << debug_step_;
    char debug_step_c[10];
    sprintf(debug_step_c, "%08d.txt", debug_step_);
    std::string debug_step_str(debug_step_c);
    grid_velocity_->Write(debug_grid_velocity_prefix_ + debug_step_str, true);
    //grid_phi_->Write(debug_grid_phi_prefix_ + debug_step_str);
    grid_pressure_->Write(debug_grid_phi_prefix_ + debug_step_str);
  }
  debug_step_++;
}  // SaveDebugDataToFile

template<common::Index N1, common::Index N2, typename T>
void FlipSimulation<N1, N2, T>::NextFrame() {
  frame_++;
}  // NextFrame

template class FlipSimulation<3,3,float>;

}  // namespace application
