#ifndef COMMON_METAGRID_H
#define COMMON_METAGRID_H

#include <cereal/archives/xml.hpp>
#include "canary/canary.h"
#include "common/primitives.h"

namespace common {

class MetaGrid {
  public:
    MetaGrid();
    MetaGrid(int x, int y, int z);
    MetaGrid(Coord global_dims);
    MetaGrid(Coord global_dims, Coord partitions, int rank);

    template<class Archive>
    void save(Archive &archive) const {
      archive(kMagicNumber);
      archive(global_dims_[0], global_dims_[1], global_dims_[2]);
      archive(partitions_[0], partitions_[1], partitions_[2]);
      archive(rank_);
      archive(kMagicNumber);
    }

    template<class Archive>
    void load(Archive &archive) {
      int magic_number = 0;
      archive(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      Coord dglobal_dims, dpartitions;
      int drank;
      archive(dglobal_dims[0], dglobal_dims[1], dglobal_dims[2]);
      archive(dpartitions[0], dpartitions[1], dpartitions[2]);
      archive(drank);
      magic_number = 0;
      archive(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      *this = MetaGrid(dglobal_dims, dpartitions, drank);
    }

  private:
    Coord global_dims_;
    Coord local_dims_;
    Coord start_;
    Coord end_;
    Coord partitions_;
    int total_partitions_;
    int rank_;
    Coord partition_id_;
    bool partitioned_;

  public:
    inline Coord global_dims() const { return global_dims_; }
    inline Coord local_dims() const { return local_dims_; }
    inline int dims(size_t i) const { return local_dims_[i]; }
    inline Coord start() const { return start_; }
    inline Coord end() const { return end_; }
    inline CoordBBox bbox() const { return CoordBBox(start_, end_); }
    inline Coord partitions() const { return partitions_; }
    inline int rank() const { return rank_; }
    inline Coord partition_id() const { return partition_id_; }
    inline bool is_partitioned() const { return partitioned_; }
    std::string string() const;

    void ComputeNeighbors(std::vector<Coord> &neighbor_ids,
                          std::vector<int> &neighbor_ranks) const;

    Coord ComputePartitionID(int rank) const;
    int ComputeRank(const Coord id) const;
    CoordBBox ComputeBBox(const Coord id) const;
    bool IsValidPartitionID(const Coord dir) const;
    bool IsValidNeighbor(int i, int j, int k) const {
      return IsValidNeighbor(Coord(i,j,k));
    }
    bool IsValidNeighbor(const Coord id) const;

    int NeighborCount() const;
};  // struct MetaGrid

}  // namespace common

#endif  // COMMON_METAGRID_H
