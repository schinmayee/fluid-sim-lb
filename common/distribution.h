#ifndef COMMON_DISTRIUBTION_H
#define COMMON_DISTRIUBTION_H

#include <cereal/archives/xml.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "canary/canary.h"
#include "common/metagrid.h"
#include "common/primitives.h"

namespace common {

class ProfileSliceLocal {
  public:
    ProfileSliceLocal();
    ProfileSliceLocal(CoordBBox local_partitions);

    void Initialize(CoordBBox local_partitions);
    void Clear();

    // Cereal archive serialization.
    template <class Archive>
    void save(Archive &archive) const {
      VLOG(1) << "Serializing dist ";
      archive(kMagicNumber);
      Coord start = local_partitions_.min();
      Coord end = local_partitions_.max();
      assert(int(count_.size()) == width_[0]*width_[1]*width_[2]);
      archive(start[0], start[1], start[2]);
      archive(end[0], end[1], end[2]);
      archive(count_);
      archive(compute_time_);
      archive(scatter_gather_time_);
      archive(kMagicNumber);
      VLOG(1) << "Deserializing dist ";
    }  // save
    template <class Archive>
    void load(Archive &archive) {
      VLOG(1) << "Deserializing dist ";
      int magic_number = 0;
      archive(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      Coord start, end;
      archive(start[0], start[1], start[2]);
      archive(end[0], end[1], end[2]);
      archive(count_);
      archive(compute_time_);
      archive(scatter_gather_time_);
      local_partitions_ = CoordBBox(start, end);
      width_ = local_partitions_.dim();
      magic_number = 0;
      archive(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      VLOG(1) << "Deserialized dist ";
    }  // load

    CoordBBox get_local_partitions() const { return local_partitions_; }

    int get_count(Coord pidx) const;
    int get_count(int i, int j, int k) const;
    void set_count(Coord pidx, int count);
    void set_count(int i, int j, int k, int count);

    void set_compute_time(float compute_time) {
      compute_time_ = compute_time;
    }
    float get_compute_time() const {
      return compute_time_;
    }

    void set_scatter_gather_time(float scatter_gather_time) {
      scatter_gather_time_ = scatter_gather_time;
    }
    float get_scatter_gather_time() const {
      return scatter_gather_time_;
    }

    std::string Print(const MetaGrid *mg) const;

  private:
    CoordBBox local_partitions_;
    Coord width_;
    std::vector<int> count_;
    float compute_time_;
    float scatter_gather_time_;

    int ComputeIndex(int i, int j, int k) const;
};  // class ProfileSliceLocal

/*
 * Profile slice for a time step.
 */
class ProfileSlice {
  public:
    ProfileSlice(const MetaGrid &mg, int total_workers,
                 int processors_per_worker);

    // Access data in this slice.

    float get_time() const { return time_; }
    void set_time(float t) { time_ = t; }

    float get_total_compute_time() const { return total_compute_time_; }
    void set_total_compute_time(float ct) { total_compute_time_ = ct; }
    float get_total_scatter_gather_time() const { return total_scatter_gather_time_; }
    void set_total_scatter_gather_time(float sgt) { total_scatter_gather_time_ = sgt; }

    int get_count(int i, int j, int k) const {
      return count_[metagrid_.ComputeRank(common::Coord(i,j,k))];
    }
    int get_count(int idx) const {
      return count_[idx];
    }
    float get_compute_time(int i, int j, int k) const {
      return compute_time_[metagrid_.ComputeRank(common::Coord(i,j,k))];
    }
    float get_compute_time(int idx) const {
      return compute_time_[idx];
    }
    void set_count(int i, int j, int k, int count) {
      int idx = metagrid_.ComputeRank(common::Coord(i,j,k));
      assert(idx >= 0);
      assert(idx < partitions_[0]*partitions_[1]*partitions_[2]);
      count_[idx] = count;
    }
    void set_compute_time(int i, int j, int k, float ct) {
      int idx = metagrid_.ComputeRank(common::Coord(i,j,k));
      assert(idx >= 0);
      assert(idx < partitions_[0]*partitions_[1]*partitions_[2]);
      compute_time_[idx] = ct;
    }
    void set_partition_stats(int i, int j, int k, int count, float ct) {
      int idx = metagrid_.ComputeRank(common::Coord(i,j,k));
      assert(idx >= 0);
      assert(idx < partitions_[0]*partitions_[1]*partitions_[2]);
      count_[idx] = count;
      compute_time_[idx] = ct;
    }

    void Gather(const std::vector<ProfileSliceLocal> &slices);

    // Compute linearized partitioning for this time.
    void ComputeLinearizedPartitionPlacement(
      canary::PartitionPlacement *placement, int offset) const;
    // Compute random partitioning, for stress testing.
    void ComputeRandomizedPartitionPlacement(
      canary::PartitionPlacement *placement, int offset) const;
    // Compute greedy partition placement.
    void ComputeGreedyPartitionPlacement(
      canary::PartitionPlacement *placement, int offset) const;

  private:
    int total_workers_;
    int processors_per_worker_;
    Coord partitions_;
    MetaGrid metagrid_;
    float time_;  // epoch -- time in simulation
    float total_compute_time_;  // compute time total
    float total_scatter_gather_time_;  // scatter/gather time total
    std::vector<int> count_;  // number of fluid cells in different regions
    std::vector<float> compute_time_;  // compute time for each partition
};  // class ProfileSlice

/*
 * FluidDistribution class stores distribution of fluid, across given chunks,
 * through time, as the simulation evolves.
 */
class FluidDistribution {
  public:
    FluidDistribution();
    FluidDistribution(const common::Coord partitions,
                      int total_workers, int processors_per_worker,
                      float current_load_factor);

    void Initialize(const common::Coord partitions,
                    int total_workers, int processors_per_worker,
                    float current_load_factor);

    common::Coord partitions() const { return partitions_; }
    int num_slices() const { return slices_.size(); }

    ProfileSlice &AddSlice(float time);
    void AddSlice(float time, const std::vector<ProfileSliceLocal> &slices);
    const ProfileSlice& GetSlice(int i) const;
    void DeleteSlice(int i);
    float GetLastTime();
    void Clear();
    void Print() const;

    const canary::PartitionHistory &GetPartitionHistory() const {
      return *history_;
    }
    void ClearPartitionHistory() {
      history_->Clear();
    }

    // Compute partitions along time, given future predictions, assuming zero
    // cost for migrating partitions.
    void ComputeZeroMigrationOraclePartitionPlacement(
      int offset, std::string partitioning);
    // Compute a partition for a fixed time/window based on estimated
    // distribution.
    void ComputeFixedWindowOraclePartitionPlacement(
      int offset, std::string partitioning, int horizon, const Coord affinity);
    // Compute one partition placement for the entire time/window based on
    // current distribution.
    void ComputeFixedWindowCurrentPartitionPlacement(
      int offset, std::string partitioning, int horizon, const Coord affinity);

    // Required archive method.
    template <class Archive>
    void serialize(Archive &ar) {
      std::cerr << "*** Serialize for FluidDistribution unimplemented ***"
                << std::endl;
    }  // serialize

  private:
    int total_workers_;
    int processors_per_worker_;
    float current_load_factor_ = 0.1;
    Coord partitions_;
    int num_partitions_;
    MetaGrid metagrid_;
    std::vector<ProfileSlice> slices_;
    canary::PartitionHistory *history_;
    // Slices available
    int steps_recorded_ = 0;
    // Steps/calls since previous partition assignment call
    int steps_since_last_reassignment_ = -1;
    // Steps assignment has been computed for (useful for fixed window)
    int steps_assigned_till_ = -1;

    // Helper methods.
    void ComputeSingleGreedyAffinePartitionPlacement(
      canary::PartitionPlacement *placement, int offset, int horizon,
      int ax, int ay, int az) const;
    void AssignSingleGreedyAffineHelper(
      const std::vector< std::vector<int> > &partition_counts,
      const std::vector<float> &partition_mean,
      int ax, int ay, int az,
      std::vector< std::vector<int> > &processor_load,
      std::vector<int> &partition_to_processor) const;

     int GetAffineCore(const Coord p, int ax, int ay, int az,
                       int num_cores) const;
     int GetAffinityCost(const Coord p, int ax, int ay, int az,
                         int num_cores, int core) const;
};  // class FluidDistribution

}  // namespace common

#endif  // COMMON_DISTRIUBTION_H
