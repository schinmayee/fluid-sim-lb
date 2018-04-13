#include <algorithm>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <utility>

#include "common/distribution.h"

#include "common/metagrid.h"
#include "common/primitives.h"
#include "common/utils.h"

namespace common {

ProfileSliceLocal::ProfileSliceLocal() {}

ProfileSliceLocal::ProfileSliceLocal(CoordBBox local_partitions) {
  Initialize(local_partitions);
}  // ProfileSliceLocal

void ProfileSliceLocal::Initialize(CoordBBox local_partitions) {
  local_partitions_ = local_partitions;
  width_ = local_partitions.dim();
  int num_partitions = width_[0] * width_[1] * width_[2];
  count_.resize(num_partitions);
  for (int i = 0; i < num_partitions; ++i) {
    count_[i] = 0;
  }
}  // Initialize

void ProfileSliceLocal::Clear() {
  compute_time_ = 0;
  scatter_gather_time_ = 0;
  width_ = local_partitions_.dim();
  int num_partitions = width_[0] * width_[1] * width_[2];
  for (int i = 0; i < num_partitions; ++i) {
    count_[i] = 0;
  }
}  // Clear

int ProfileSliceLocal::get_count(Coord pidx) const {
  get_count(pidx[0], pidx[1], pidx[2]);
}  // get_count

int ProfileSliceLocal::get_count(int i, int j, int k) const {
  return count_[ComputeIndex(i, j, k)];
}  // get_count

void ProfileSliceLocal::set_count(Coord pidx, int count) {
  set_count(pidx[0], pidx[1], pidx[2], count);
}  // set_count

void ProfileSliceLocal::set_count(int i, int j, int k, int count) {
  count_[ComputeIndex(i, j, k)] = count;
}  // set_count

int ProfileSliceLocal::ComputeIndex(int i, int j, int k) const {
  const Coord start = local_partitions_.min();
  int di = i - start[0];
  int dj = j - start[1];
  int dk = k - start[2];
  CHECK(di >= 0 && di < width_[0]);
  CHECK(dj >= 0 && dj < width_[1]);
  CHECK(dk >= 0 && dk < width_[2]);
  int idx = (dk * width_[1] + dj) * width_[0] + di;
  CHECK(idx >= 0);
  CHECK(idx < int(count_.size()));
  return idx;
}  // ComputeIndex

std::string ProfileSliceLocal::Print(const MetaGrid *mg) const {
  std::stringstream ss;
  Coord start = local_partitions_.min();
  Coord end = local_partitions_.max();
  ss << "Partitions " << ToString(local_partitions_) << std::endl;
  ss << "* Compute time : " << compute_time_ << std::endl;
  ss << "* Scatter-gather time : " << scatter_gather_time_ << std::endl;
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        Coord pid(i,j,k);
        ss << "LocalCount for " << ToString(pid);
        ss << " , count " << get_count(pid);
        if (mg) {
          ss << ", rank " << mg->ComputeRank(pid);
        }
        ss << std::endl;
      }  // for k
    }  // for j
  }  // for i
  return ss.str();
}  // Print

ProfileSlice::ProfileSlice(
    const MetaGrid &mg, int total_workers, int processors_per_worker)
  : total_workers_(total_workers),
    processors_per_worker_(processors_per_worker),
    partitions_(mg.partitions()),
    metagrid_(mg), time_(0),
    total_compute_time_(0), total_scatter_gather_time_(0) {
  int num_partitions = partitions_[0]*partitions_[1]*partitions_[2];
  count_.resize(num_partitions);
  compute_time_.resize(num_partitions);
  for (int i = 0; i < num_partitions; ++i) {
    count_[i] = 0;
    compute_time_[i] = 0;
  }  // for i < num_partitions_
}  // ProfileSlice

void ProfileSlice::Gather(const std::vector<ProfileSliceLocal> &slices) {
  const Coord partitions = metagrid_.partitions();
  int num_partitions = partitions[0]*partitions[1]*partitions[2];
  total_compute_time_ = 0;
  total_scatter_gather_time_ = 0;
  for (const auto &s : slices) {
    CoordBBox local_partitions = s.get_local_partitions();
    Coord start = local_partitions.min();
    Coord end   = local_partitions.max();
    const float compute_time = s.get_compute_time();
    int total_count = 0;
    for (int i = start[0]; i <= end[0]; ++i) {
      for (int j = start[1]; j <= end[1]; ++j) {
        for (int k = start[2]; k <= end[2]; ++k)  {
          total_count += (s.get_count(i, j, k) + 1);
        }  // for k
      }  // for j
    }  // for i
    float inv_total_count = 1.0/float(total_count);
    for (int i = start[0]; i <= end[0]; ++i) {
      for (int j = start[1]; j <= end[1]; ++j) {
        for (int k = start[2]; k <= end[2]; ++k)  {
          int count = s.get_count(i, j, k);
          float ct = compute_time * inv_total_count * float(count+1);
          set_partition_stats(i, j, k, count, ct);
        }  // for k
      }  // for j
    }  // for i
    total_compute_time_ += compute_time;
    total_scatter_gather_time_ += s.get_scatter_gather_time();
  }  // for s
}  // Gather

void ProfileSlice::ComputeLinearizedPartitionPlacement(
  canary::PartitionPlacement *placement, int offset) const {
  CHECK(partitions_ == metagrid_.partitions());
  int num_partitions = partitions_[0]*partitions_[1]*partitions_[2];
  int linearized_count[num_partitions];
  int total_count = 0;
  // Linearize the array of partition counts. Also compute total work.
  for (int i = 0; i < partitions_[0]; ++i) {
    for (int j = 0; j < partitions_[1]; ++j) {
      for (int k = 0; k < partitions_[2]; ++k) {
        int idx = metagrid_.ComputeRank(common::Coord(i,j,k));
        linearized_count[idx] = count_[idx];
        total_count += count_[idx];
      }  // for k
    }  // for j
  }  // for i
  // Compute work per processor.
  int total_workers = total_workers_ - offset;
  int count_per_processor =
    ceil(float(total_count)/float(total_workers*processors_per_worker_));
  int worker_processors_assigned = 0;
  int work_on_processor = 0;
  int next_worker = offset;
  int next_processor = 0;
  // Assign partitions to workers, assigning every
  // count_per_processor*processors_per_worker_ to a worker.
  for (int i = 0; i < num_partitions; ++i) {
    if (work_on_processor >= count_per_processor ||
        ((work_on_processor > 0.75*count_per_processor) &&
         (work_on_processor + linearized_count[i] > 1.5*count_per_processor))) {
      worker_processors_assigned++;
      work_on_processor = 0;
    }  // if processor_count ...
    if (worker_processors_assigned == processors_per_worker_) {
      worker_processors_assigned = 0;
      if (next_worker < total_workers-1) {
        next_worker++;
      }  // cap next_worker at total_workers-1
    }  // if worker_processors_assigned == processors_per_worker_
    CHECK(next_worker >= 0);
    CHECK(next_worker < total_workers);
    work_on_processor += linearized_count[i];
    VLOG(1) << "Placing partition " << i << " on worker "
              << next_worker;
    placement->PlacePartitionOnWorker(i, next_worker+offset);  // start at 1
  }  // for i < num_partitions
  VLOG(1) << "Partition placement: ";
  for (int i = 0; i < partitions_[0]; ++i) {
    for (int j = 0; j < partitions_[1]; ++j) {
      for (int k = 0; k < partitions_[2]; ++k) {
        int idx = metagrid_.ComputeRank(common::Coord(i,j,k));
        VLOG(1) << "* " << common::ToString(common::Coord(i,j,k)) << " : " <<
                     placement->GetPartitionPlacement(idx);
      }  // for k
    }  // for j
  }  // for i
}  // ComputeLinearizedPartitionPlacement

void ProfileSlice::ComputeRandomizedPartitionPlacement(
  canary::PartitionPlacement *placement, int offset) const {
  CHECK(partitions_ == metagrid_.partitions());
  int total_workers = total_workers_ - offset;
  int num_partitions = partitions_[0]*partitions_[1]*partitions_[2];
  int partitions_per_worker = ceil(float(num_partitions)/float(total_workers));
  // Create a random partition assignment.
  std::vector<int> assignment(num_partitions);
  int assigned_to_worker = 0;
  int worker_id = 0;
  for (int i = 0; i < num_partitions; ++i) {
    if (assigned_to_worker < partitions_per_worker) {
      assignment[i] = worker_id;
      assigned_to_worker += 1;
    } else {
      worker_id++;
      assignment[i] = worker_id;
      assigned_to_worker = 1;
    }
  }
  CHECK(worker_id < total_workers);
  std::random_shuffle(assignment.begin(), assignment.end());
  // Assign partitions to workers.
  for (int i = 0; i < num_partitions; ++i) {
    placement->PlacePartitionOnWorker(i, assignment[i]+offset);  // start at offset
  }  // for i < num_partitions
  VLOG(1) << "Partition placement: ";
  for (int i = 0; i < partitions_[0]; ++i) {
    for (int j = 0; j < partitions_[1]; ++j) {
      for (int k = 0; k < partitions_[2]; ++k) {
        int idx = metagrid_.ComputeRank(common::Coord(i,j,k));
        VLOG(1) << "* " << common::ToString(common::Coord(i,j,k)) << " : " <<
                     placement->GetPartitionPlacement(idx);
      }  // for k
    }  // for j
  }  // for i
}  // ComputeRandomizedPartitionPlacement

void ProfileSlice::ComputeGreedyPartitionPlacement(
  canary::PartitionPlacement *placement, int offset) const {
  int total_workers = total_workers_ - offset;
  int total_processors = total_workers * processors_per_worker_;
  int num_partitions = partitions_[0]*partitions_[1]*partitions_[2];
  // Store assignment and other quantities to compute assignment.
  std::vector<int> assignment(num_partitions);
  // Store negative work -- pick the largest element for assigning.
  // Each element stores assignment and processor id.
  std::vector<std::pair<float, int>> processor_work(total_processors);
  for (int i = 0; i < total_processors; ++i) {
    processor_work[i].first = 0;
    processor_work[i].second = i;
  }
  // Make assignment into heap so that it is easy to pull out the processor with
  // least amount of work.
  std::make_heap(processor_work.begin(), processor_work.end());
  // Sort based on compute time.
  auto sorted_indexes = SortIndexes(compute_time_);
  for (auto id : sorted_indexes) {
    float work_to_assign = compute_time_[id];
    // Pick the processor with least amount of work.
    std::pop_heap(processor_work.begin(), processor_work.end());
    // Update cost/assignment. Use negative work.
    // Also update assignment for the partition.
    assignment[id] = processor_work.back().second/processors_per_worker_;
    LOG(INFO) << "Assigning " << id << " : " << work_to_assign << " to "
              << assignment[id] << " which had previously "
              << processor_work.back().first;
    processor_work.back().first -= work_to_assign;
    // Add this back to the heap.
    std::push_heap(processor_work.begin(), processor_work.end());
  }  // for id
  // Assign partitions to workers.
  for (int i = 0; i < num_partitions; ++i) {
    placement->PlacePartitionOnWorker(i, assignment[i]+offset);  // start at offset
  }  // for i < num_partitions
  VLOG(1) << "Partition placement: ";
  for (int i = 0; i < partitions_[0]; ++i) {
    for (int j = 0; j < partitions_[1]; ++j) {
      for (int k = 0; k < partitions_[2]; ++k) {
        int idx = metagrid_.ComputeRank(common::Coord(i,j,k));
        VLOG(1) << "* " << common::ToString(common::Coord(i,j,k)) << " : " <<
                     placement->GetPartitionPlacement(idx);
      }  // for k
    }  // for j
  }  // for i
}  // ComputeGreedyPartitionPlacementUsingComputeTime

FluidDistribution::FluidDistribution() : history_(nullptr) {}

FluidDistribution::FluidDistribution(
    const Coord partitions,
    int total_workers, int processors_per_worker,
    float current_load_factor)
  : partitions_(partitions), metagrid_(partitions, partitions, 0),
    total_workers_(total_workers),
    processors_per_worker_(processors_per_worker),
    current_load_factor_(current_load_factor) {
  num_partitions_ = partitions[0]*partitions[1]*partitions[2];
  history_ = new canary::PartitionHistory(num_partitions_);
  srand(0);
}  // FluidDistribution

void FluidDistribution::Initialize(
    const Coord partitions,
    int total_workers, int processors_per_worker,
    float current_load_factor) {
  partitions_ = partitions;
  metagrid_ = MetaGrid(partitions, partitions, 0);
  total_workers_ = total_workers;
  processors_per_worker_ = processors_per_worker;
  current_load_factor_ = current_load_factor;
  num_partitions_ = partitions[0]*partitions[1]*partitions[2];
  history_ = new canary::PartitionHistory(num_partitions_);
  srand(0);
}  // Initialize

ProfileSlice& FluidDistribution::AddSlice(float time) {
  CHECK(total_workers_ > 0);
  CHECK(processors_per_worker_ > 0);
  slices_.emplace_back(metagrid_, total_workers_, processors_per_worker_);
  ProfileSlice &slice = slices_.back();
  slice.set_time(time);
  steps_recorded_++;
  return slices_.back();
}  // AddSlice

void FluidDistribution::AddSlice(float time, const std::vector<ProfileSliceLocal> &slices) {
  ProfileSlice &s = AddSlice(time);
  s.Gather(slices);
}  // AddSlice

const ProfileSlice& FluidDistribution::GetSlice(int i) const {
  CHECK(i >= 0);
  CHECK(i < int(slices_.size()));
  return slices_[i];
}  // GetSlice

void FluidDistribution::DeleteSlice(int i) {
  CHECK(i >= 0);
  CHECK(i < int(slices_.size()));
  slices_.erase(slices_.begin() + i);
}  // DeleteSlice

float FluidDistribution::GetLastTime() {
  CHECK(slices_.size() > 0);
  return slices_.back().get_time();
}  // GetLastTime

void FluidDistribution::Clear() {
  slices_.clear();
}  // Clear

void FluidDistribution::Print() const {
  int num_slices = slices_.size();
  if (num_slices == 0) {
    return;
  }
  const ProfileSlice &slice = slices_[num_slices-1];
  LOG(INFO) << "Time " << slice.get_time() << ", " <<
               "ComputeTime " << slice.get_total_compute_time() << ", " <<
               "ScatterGatherTime " << slice.get_total_scatter_gather_time();
  for (int i = 0; i < partitions_[0]; ++i) {
    for (int j = 0; j < partitions_[1]; ++j) {
      for (int k = 0; k < partitions_[2]; ++k) {
        LOG(INFO) << "* Partition " << i << " " << j << " " << k << ", " <<
                     "Count " << slice.get_count(i, j, k) << ", " <<
                     "ComputeTime " << slice.get_compute_time(i, j, k);
      }  // for k
    }  // for j
  }  // for i
  // Print count per worker using most recently computed placement.
  const canary::PartitionPlacement &previous_placement =
    history_->PreviousPlacement();
  int num_partitions = previous_placement.GetNumPartitions();
  CHECK_EQ(num_partitions, num_partitions_);
  std::vector<int> worker_count(total_workers_, 0);
  int total_count = 0;
  for (int i = 0; i < num_partitions; ++i) {
    int wid = previous_placement.GetPartitionPlacement(i);
    int count = slice.get_count(i);
    worker_count[wid] += count;
    total_count += count;
  }  // i < num_partitions
  int total_workers = (worker_count[0] == 0) ? total_workers_ - 1 : total_workers_;
  int average = total_count/total_workers_ + 1;
  LOG(INFO) << "Counts per worker:";
  for (int w = 0; w < total_workers_; ++w) {
    LOG(INFO) << "* Worker " << w <<
                 " Count " << worker_count[w];
  }
  int maximum = *std::max_element(
    std::begin(worker_count), std::end(worker_count));
  LOG(INFO) << "Average worker count = " << average <<
              " and maximum = " << maximum;
}  // Print

// Greedy in time, ignores migration cost.
void FluidDistribution::ComputeZeroMigrationOraclePartitionPlacement(
  int offset, std::string partitioning) {
  if (slices_.size() == 0) {
    return;
  }
  const canary::PartitionPlacement &previous_placement = history_->PreviousPlacement();
  int ns = num_slices();
  for (int i = 0; i < ns; ++i) {
    canary::PartitionPlacement placement(previous_placement.GetNumPartitions());
    const ProfileSlice &slice = GetSlice(i);
    VLOG(1) << "Time for new placement : " << slice.get_time();
    if (partitioning == "linear") {
      slice.ComputeLinearizedPartitionPlacement(&placement, offset);
    } else if (partitioning == "random") {
      slice.ComputeRandomizedPartitionPlacement(&placement, offset);
    } else if (partitioning == "greedy") {
      slice.ComputeGreedyPartitionPlacement(&placement, offset);
    } else {
      CHECK(false) << "Invalid partitioning " << partitioning;
    }
    if (previous_placement.IsInvalid() || !placement.IsSameAs(previous_placement)) {
      history_->AddPartitionPlacement(slice.get_time(), placement);
    }  // if !placement.IsSameAs(previous_placement)
  }  // for i < ns
  // Add last time regardless.
  history_->SetLastTime(GetSlice(ns-1).get_time());
}  // ComputeZeroMigrationOraclePartitionPlacement

int FluidDistribution::GetAffineCore(
  const Coord p, int ax, int ay, int az, int num_cores)  const {
  const Coord partitions = metagrid_.partitions();
  int lin_per_core = ax*ay*az/num_cores;
  const int nx = partitions[0]/ax;
  const int ny = partitions[1]/ay;
  const int nz = partitions[2]/az;
  int c0 = p[0]/nx;
  int c1 = p[1]/ny;
  int c2 = p[2]/nz;
  int c_lin = c0 + c1*ax + c2*ax*ay;
  VLOG(1) << "LinID/Affinity " << ToString(p) << " is " << c_lin
          << "/" << (c_lin/lin_per_core);
  return c_lin/lin_per_core;
}  // GetAffineCore

int FluidDistribution::GetAffinityCost(
  const Coord p, int ax, int ay, int az, int num_cores, int core) const {
  if (GetAffineCore(p, ax, ay, az, num_cores) == core) {
    return 0;
  } else {
    return 1;
  }
}  // GetAffinityCost

void FluidDistribution::AssignSingleGreedyAffineHelper(
  const std::vector< std::vector<int> > &partition_counts,
  const std::vector<float> &partition_mean,
  int ax, int ay, int az,
  std::vector< std::vector<int> > &processor_load,
  std::vector<int> &partition_to_processor) const {
  int num_partitions = partition_counts.size();
  int num_cores = processor_load.size();
  int horizon = partition_counts[0].size();
  // Group partitions into almost empty and not empty.
  // Want to distribute almost empty partitions evenly, rather than all of them
  // going to a single node, since there is a fixed overhead associated with
  // each.
  float max_partition_load = *std::max_element(
    std::begin(partition_mean), std::end(partition_mean));
  int thresh = 0.04*max_partition_load;
  std::vector<int> non_empty_idxs;
  std::vector<int> empty_idxs;
  for (int p = 0; p < num_partitions; ++p) {
    if (partition_mean[p] <= thresh) {
      empty_idxs.push_back(p);
    } else {
      non_empty_idxs.push_back(p);
    }
  }  // for p
  // Get sorted indices for non-empty partitions.
  std::sort(non_empty_idxs.begin(), non_empty_idxs.end(), 
            [this, &partition_mean](size_t i1, size_t i2) {
              if (partition_mean[i1] != partition_mean[i2]) {
                return partition_mean[i1] > partition_mean[i2]; 
              }
              Coord c1 = metagrid_.ComputePartitionID(i1);
              Coord c2 = metagrid_.ComputePartitionID(i2);

              if (c1.x() != c2.x()) {
                return c1.x() > c2.x();
              }

              if (c1.y() != c2.y()) {
                return c1.y() > c2.y();
              }

              return c1.z() > c2.z();
            });
  // Get average load at each step, this is the best case.
  float inv_cores = 1.0/float(num_cores);
  std::vector<int> ave_load(horizon, 0);
  std::vector<float> max_load(horizon, 0);
  for (int t = 0; t < horizon; ++t) {
    for (int p = 0; p < num_partitions; ++p) {
      ave_load[t] += partition_counts[p][t];
    }
    float ave_load_t = inv_cores*ave_load[t];
    ave_load[t] = ave_load_t;
  }
  // Max for each time step.
  //std::vector<int> max_counts(horizon, 0);
  const float factor = 0.5;
  const int affinity_cost = 10;
  // Go over non-empty partitions and assign them.
  std::vector<int> processor_num_partitions(num_cores, 0);
  for (const auto p : non_empty_idxs) {
    const std::vector<int> &p_counts = partition_counts[p];
    Coord pcoord = metagrid_.ComputePartitionID(p);
    // First check the cost of assigning this partition to the affine core.
    // If it does not exceed the ideal cost, assign it to the  affine core.
    int affine_proc = GetAffineCore(pcoord, ax, ay, az, num_cores);
    int cost_to_affine = 0;
    int best_cost = 0;
    for (int t = 0; t < horizon; ++t) {
      int load_t = processor_load[affine_proc][t] + p_counts[t];
      cost_to_affine += std::max(load_t, ave_load[t]);
      best_cost += ave_load[t];
    }  // for t < horizon
    int proc_id = -1;
    if (cost_to_affine < best_cost + horizon * affinity_cost) {
      proc_id = affine_proc;
    } else {
      // Compute the cost of assgining this partition to each possible core.
      std::vector<float> assignment_cost(num_cores, 0);
      for (int c = 0; c < num_cores; ++c) {
        for (int t = 0; t < horizon; ++t) {
          // Max factor
          assignment_cost[c] += (std::max(max_load[t], float(processor_load[c][t]))
                                 + p_counts[t]);
          // Current load factor -- don't want to put too many with the same
          // trend on the same node. If there are multiple candidates that do
          // not exeed max at each time in the horizon, pick the candiate with
          // the least amount of work.
          assignment_cost[c] += factor*(processor_load[c][t] + p_counts[t]);
        }  // for t < horizon
      }  // for c
      auto min_elem = std::min_element(
        std::begin(assignment_cost), std::end(assignment_cost));
      proc_id = std::distance(std::begin(assignment_cost), min_elem);
    }  // if affine else
    // Finalize the assignment -- update processor load and max for each t.
    partition_to_processor[p] = proc_id;
    processor_num_partitions[proc_id]++;
    LOG(INFO) << "Assigning partition " << p << " " <<
                  ToString(metagrid_.ComputePartitionID(p)) <<
                 " load " << partition_mean[p] << " to " << proc_id;
    //processor_mean_load[proc_id] += partition_mean[p];
    for (int t = 0; t < horizon; ++t) {
      processor_load[proc_id][t] += partition_counts[p][t];
      max_load[t] = std::max(max_load[t], float(processor_load[proc_id][t]));
    }  // for t
  }  // for p
  // Go over partitions that are almost empty and generate an initial assignment
  // for them.
  int num_empty = empty_idxs.size();
  std::vector<int> processor_num_partitions_temp = processor_num_partitions;
  std::vector<std::vector<int>> processor_empty_partitions(num_cores);
  for (int i = 0; i < num_empty; ++i) {
    int partition_id = empty_idxs[i];
    Coord pcoord = metagrid_.ComputePartitionID(partition_id);
    int proc_id = GetAffineCore(pcoord, ax, ay, az, num_cores);
    processor_num_partitions_temp[proc_id]++;
    processor_empty_partitions[proc_id].push_back(partition_id);
  }
  // Now reassign excess to equalize total partitions.
  int max_partitions_per_proc = num_partitions/num_cores + 4;
  LOG(INFO) << "Maximum partitions per proc should be " << max_partitions_per_proc;
  for (int c = 0; c < num_cores; ++c) {
    int reassign_num = processor_num_partitions_temp[c] - 
                       std::max(max_partitions_per_proc, processor_num_partitions[c]);
    if (reassign_num > 0) {
      for (int i = 0; i < reassign_num; ++i) {
        auto min_elem = std::min_element(
          std::begin(processor_num_partitions_temp),
          std::end(processor_num_partitions_temp));
        int to_c = std::distance(std::begin(processor_num_partitions_temp), min_elem);
        processor_num_partitions_temp[c]--;
        processor_num_partitions_temp[to_c]++;
        processor_num_partitions[to_c]++;
        processor_empty_partitions[to_c].push_back(processor_empty_partitions[c].back());
        processor_empty_partitions[c].pop_back();
      }
    }
  }
  // Finalize this assignment.
  for (int c = 0; c < num_cores; ++c) {
    for (auto p : processor_empty_partitions[c]) {
      partition_to_processor[p] = c;
      for (int t = 0; t < horizon; ++t) {
        processor_load[c][t] += partition_counts[p][t];
      }  // for t
      LOG(INFO) << "Assigning partition " << p << " " <<
                    ToString(metagrid_.ComputePartitionID(p)) <<
                   " load " << partition_mean[p] << " to " << c;
    }
  }
}  // AssignSingleGreedyAffineHelper

void FluidDistribution::ComputeSingleGreedyAffinePartitionPlacement(
  canary::PartitionPlacement *placement, int offset, int horizon,
  int ax, int ay, int az) const {
  // Compute min, max, and variation, max variation and mean load for each
  // partition.
  std::vector<float> mean_counts(num_partitions_, 0);
  std::vector< std::vector<int> > counts(num_partitions_);
  for (int p = 0; p < num_partitions_; ++p) {
    std::vector<int> cnts(horizon, 0);
    for (int s = 0; s < horizon; ++s) {
      int c = slices_[s].get_count(p);
      cnts[s] = c;
    }  // for s
    counts[p] = cnts;
    int total_count = std::accumulate(std::begin(cnts), std::end(cnts), 0);
    mean_counts[p] = float(total_count)/float(horizon);
  }  // for p
  // Total number of processors.
  int num_workers = total_workers_ - offset;
  int num_processors = num_workers * processors_per_worker_;
  CHECK_EQ((ax*ay*az)%num_workers, 0);
  // Assign all partitions using a simple greedy algorithm.
  std::vector<int> partition_to_processor(num_partitions_);
  std::vector<std::vector<int>> processor_load(
    num_processors, std::vector<int>(horizon, 0));
  AssignSingleGreedyAffineHelper(
    counts, mean_counts, ax, ay, az, processor_load, partition_to_processor);
  // Now use computed assignment to assign partitions to workers.
  // Round-robin across workers so that there is more even placement.
  
  std::unordered_map<int, std::vector<float>> proc_to_mean_counts;
  for (int p = 0; p < num_partitions_; ++p) {
    int proc = partition_to_processor[p];
    proc_to_mean_counts[proc].push_back(mean_counts[p]);
    int wid = proc % num_workers;
    placement->PlacePartitionOnWorker(p, wid+offset);
  }

  for (int wid = 0; wid < num_workers; ++wid) {
    std::stringstream log;
    log <<  "wid: " << wid << "\n";
    for (int proc_id = wid; proc_id < num_processors; proc_id += num_workers) {
      auto it = proc_to_mean_counts.find(proc_id);
      log << "proc: " << proc_id;
      if (it != proc_to_mean_counts.end()) {
        std::sort(it->second.begin(), it->second.end());
        float sum = 0;
        for (auto mc : it->second) {
          sum += mc;
        }
        log << "n: " << it->second.size() << " s: " << sum << " |";
        for (auto mc : it->second) {
          log << " " << mc;
        }
      }
      log << "\n";
    }
    LOG(INFO) << log.str();
  }
}  // ComputeSingleOraclePartitionPlacement

void FluidDistribution::ComputeFixedWindowOraclePartitionPlacement(
    int offset, std::string partitioning, int horizon, const Coord affinity) {
  if (steps_since_last_reassignment_ == -1) {
    steps_since_last_reassignment_ = 1;
  }
  if (steps_since_last_reassignment_ >= horizon) {
    int ns = num_slices();
    CHECK(ns >= horizon) << "Not enough slices to compute new partitioning, " <<
                            "only " << ns << " avaiable, need " << horizon;
    const canary::PartitionPlacement &previous_placement =
      history_->PreviousPlacement();
    canary::PartitionPlacement placement(previous_placement.GetNumPartitions());
    LOG(INFO) << "Computing new partition assignment ...";
    ComputeSingleGreedyAffinePartitionPlacement(
      &placement, offset, horizon,
      affinity[0], affinity[1], affinity[2]);
    steps_since_last_reassignment_ = 1;
    steps_assigned_till_ += horizon;
    // Add this single placement to history.
    history_->Clear();
    // This partitioning can be activated starting time slices_[0].get_time()
    history_->AddPartitionPlacement(slices_[0].get_time(), placement);
    // The recorded placements are for time till slices_[horizon-1].get_time()
    history_->SetLastTime(slices_[horizon-1].get_time());
    // Remove slices used for next window.
    slices_.erase(slices_.begin(), slices_.begin() + horizon);
  } else {
    LOG(INFO) << "Not computing new partition assignment ...";
    steps_since_last_reassignment_++;
  }
}  // ComputeFixedOraclePartitionPlacement

void FluidDistribution::ComputeFixedWindowCurrentPartitionPlacement(
    int offset, std::string partitioning, int horizon, const Coord affinity) {
  if (slices_.size() == 0) {
    LOG(INFO) << "No slices to compute new partitioning";
    return;
  }
  const canary::PartitionPlacement &previous_placement = history_->PreviousPlacement();
  int ns = num_slices();
  const ProfileSlice &slice = GetSlice(ns-1);
  if (steps_since_last_reassignment_ == -1 ||
      steps_since_last_reassignment_ >= horizon) {
    LOG(INFO) << "Computing new partition assignment ...";
    steps_since_last_reassignment_ = 1;
    canary::PartitionPlacement placement(previous_placement.GetNumPartitions());
    if (partitioning == "linear") {
      slice.ComputeLinearizedPartitionPlacement(&placement, offset);
    } else if (partitioning == "random") {
      slice.ComputeRandomizedPartitionPlacement(&placement, offset);
    } else if (partitioning == "greedy") {
      ComputeSingleGreedyAffinePartitionPlacement(
        &placement, offset, 1,
        affinity[0], affinity[1], affinity[2]);
      //slice.ComputeGreedyPartitionPlacement(&placement, offset);
    } else {
      CHECK(false) << "Invalid partitioning " << partitioning;
    }
    history_->Clear();
    history_->AddPartitionPlacement(slice.get_time(), placement);
  } else {
    LOG(INFO) << "Not computing new partition assignment ...";
    steps_since_last_reassignment_++;
  }
  history_->SetLastTime(slice.get_time());
}  // ComputeFixedCurrentPartitionPlacement

}  // namespace common
