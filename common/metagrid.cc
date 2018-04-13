#include "common/metagrid.h"
#include "common/primitives.h"


namespace common {

MetaGrid::MetaGrid()
  : global_dims_(0), local_dims_(0),
    start_(0), end_(0), partitions_(1),
    total_partitions_(partitions_[0]*partitions_[1]*partitions_[2]),
    rank_(0), partition_id_(0), partitioned_(false) {}

MetaGrid::MetaGrid(int x, int y, int z)
  : global_dims_(x,y,z), local_dims_(x,y,z),
    start_(0), end_(0), partitions_(1),
    total_partitions_(partitions_[0]*partitions_[1]*partitions_[2]),
    rank_(0), partition_id_(0), partitioned_(false) {}

MetaGrid::MetaGrid(Coord global_dims)
  : global_dims_(global_dims), local_dims_(global_dims),
    start_(0), end_(0), partitions_(1),
    total_partitions_(partitions_[0]*partitions_[1]*partitions_[2]),
    rank_(0), partition_id_(0), partitioned_(false) {}

MetaGrid::MetaGrid(Coord global_dims, Coord partitions, int rank)
  : global_dims_(global_dims), partitions_(partitions),
    total_partitions_(partitions_[0]*partitions_[1]*partitions_[2]),
    rank_(rank), partitioned_(true) {
  partition_id_ = ComputePartitionID(rank);
  for (size_t d = 0; d < 3; ++d) {
    local_dims_[d] = global_dims[d]/partitions_[d];
    start_[d] = local_dims_[d] * partition_id_[d];
    end_[d] = start_[d] + local_dims_[d] - 1;
  }  // for d < 3
}  // MetaGrid

std::string MetaGrid::string() const {
  std::string result = "";
  result += ("Global dimensions: " + ToString(global_dims()) + "\n");
  result += ("Rank: " + std::to_string(rank()) + "\n");
  result += ("Partition id: " + ToString(partition_id()) + "\n");
  result += ("Local dimensions: " + ToString(local_dims()) + "\n");
  result += ("Start: " + ToString(start()) + "\n");
  result += ("End: " + ToString(end()) + "\n");
  return result;
}  // string

void MetaGrid::ComputeNeighbors(
  std::vector<Coord> &neighbor_ids, std::vector<int> &neighbor_ranks) const {
  neighbor_ids.clear();
  neighbor_ranks.clear();
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        if (i == 0 && j == 0 && k == 0) {
          continue;
        }
        Coord nid = partition_id_ + Coord(i,j,k);
        if (nid[0] < 0 || nid[0] >= partitions_[0] ||
            nid[1] < 0 || nid[1] >= partitions_[1] ||
            nid[2] < 0 || nid[2] >= partitions_[2]) {
          continue;
        }
        neighbor_ids.push_back(nid);
        neighbor_ranks.push_back(ComputeRank(nid));
      }  // for k
    }  // for j
  }  // for i
  CHECK_EQ(neighbor_ids.size(), neighbor_ranks.size());
}  // ComputeNeighbors

Coord MetaGrid::ComputePartitionID(int rank) const {
  CHECK(partitions_[0] != 0);
  CHECK(partitions_[1] != 0);
  CHECK(partitions_[2] != 0);
  Coord id;
  id[0] = rank % partitions_[0];
  id[1] = ((rank - id[0])/partitions_[0]) % partitions_[1];
  id[2] = (rank - id[0] - id[1]*partitions_[0])/(partitions_[0]*partitions_[1]);
  CHECK_EQ(ComputeRank(id), rank) << "Error in rank/partition id computation";
  return id;
}  // ComputePartitionID

int MetaGrid::ComputeRank(const Coord id) const {
  int rank = (id[0] + partitions_[0]*(id[1] + partitions_[1]*id[2]));
  return rank;
}  // ComputeRank

CoordBBox MetaGrid::ComputeBBox(const Coord id) const {
  Coord partition_start, partition_end;
  for (size_t d = 0; d < 3; ++d) {
    partition_start[d] = id[d] * local_dims_[d];
    partition_end[d] = partition_start[d] + local_dims_[d] - 1;
  }  // for d < 3
  return CoordBBox(partition_start, partition_end);
}  // ComputeBounds

bool MetaGrid::IsValidPartitionID(const Coord id) const {
  return !(id[0] < 0 || id[0] >= partitions_[0] ||
           id[1] < 0 || id[1] >= partitions_[1] ||
           id[2] < 0 || id[2] >= partitions_[2]);
}  // IsValidPartitionID

bool MetaGrid::IsValidNeighbor(const Coord dir) const {
  Coord neighbor_id = partition_id_ + dir;
  return IsValidPartitionID(neighbor_id);
}  // IsValidNeighbor

int MetaGrid::NeighborCount() const {
  int count = 0;
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        if (i == 0 && j == 0 && k == 0) {
          continue;
        }
        Coord neighbor_id = partition_id_ + Coord(i,j,k);
        if (IsValidPartitionID(neighbor_id)) {
          count++;
        }
      }  // for k
    }  // for j
  }  // for i
  return count;
}  // NeighborCount

}  // namespace common
