#include "canary/canary.h"
#include "common/vector.h"
#include "common/metagrid.h"
#include "common/scalar_grid.h"

namespace common {

template class Vector<float>;
template class Vector<double>;

template<typename T> template<class ScalarGridVIdx>
void Vector<T>::SendGhostData(
  canary::CanaryTaskContext *task_context, int ghost_width,
  const ScalarGridVIdx &idx) const {
  int count = 0;
  const MetaGrid metagrid = idx.metagrid();
  // Get neighbor rank and serialize data for each neighbor.
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        if (i == 0 && j == 0 && k == 0) {
          continue;
        }
        if (metagrid.IsValidNeighbor(i,j,k)) {
          struct evbuffer *buffer = evbuffer_new();
          canary::CanaryOutputArchive archive(buffer);
          int rank = SerializeGhost(
            i, j, k, ghost_width, idx, archive);
          assert(rank >= 0);
          task_context->Scatter(rank, canary::RawEvbuffer(buffer));
          count++;
        }  // if valid neighbor
      }  // for k
    }  // for j
  }  // for i
  assert(count == metagrid.NeighborCount());
}  // SendGhostData

template<typename T> template<class ScalarGridVIdx>
int Vector<T>::ReceiveGhostData(
  canary::CanaryTaskContext *task_context, int ghost_width,
  const ScalarGridVIdx &idx) {
  const MetaGrid metagrid = idx.metagrid();
  // Go over data rceived in receive buffer.
  EXPECT_GATHER_SIZE(metagrid.NeighborCount());
  auto recv_buffer = task_context->Gather<canary::RawEvbuffer>();
  for (auto &raw_buffer : recv_buffer) {
    canary::CanaryInputArchive archive(raw_buffer.buffer);
    DeserializeGhost(ghost_width, idx, archive);
    evbuffer_free(raw_buffer.buffer);
  }  // for raw_buffer
  return 0;
}  // ReceiveGhostData

template<typename T> template<class ScalarGridVIdx>
int Vector<T>::SerializeGhost(
  int di, int dj, int dk, int ghost_width, const ScalarGridVIdx &idx,
  canary::CanaryOutputArchive &archive) const {
  return SerializeGhost(Coord(di,dj,dk), ghost_width, idx, archive);
}  // SerializeGhost

template<typename T> template<class ScalarGridVIdx>
int Vector<T>::SerializeGhost(
  Coord neighbor_dir, int ghost_width, const ScalarGridVIdx &idx,
  canary::CanaryOutputArchive &archive) const {
  const MetaGrid metagrid = idx.metagrid();
  const Coord pid = metagrid.partition_id();
  const Coord neighbor_id = pid + neighbor_dir;
  if (!metagrid.IsValidPartitionID(neighbor_id)) {
    VLOG(1) << "No neighbor in " << ToString(neighbor_dir) << " for "
              << ToString(pid);
    return -1;
  }
  CoordBBox bbox = metagrid.ComputeBBox(neighbor_id);
  VLOG(1) << "Original box for neighbor " << ToString(neighbor_dir) <<
               " of " << ToString(metagrid.partition_id()) << " is " <<
               ToString(bbox);
  bbox.expand(ghost_width);
  bbox.intersect(metagrid.bbox());
  VLOG(1) << "Box for neighbor " << ToString(neighbor_dir) << " of " <<
               ToString(metagrid.partition_id()) << " is " <<
               ToString(bbox); 
  SerializeDense(bbox, idx, archive);
  int rank = metagrid.ComputeRank(neighbor_id);
  return rank;
}  // SerializeGhost

template<typename T> template<class ScalarGridVIdx>
void Vector<T>::SerializeDense(
  const CoordBBox box, const ScalarGridVIdx &idx,
  canary::CanaryOutputArchive &archive) const {
  typedef typename ScalarGridVIdx::ValueType VIdx;
  const Coord start = box.min();
  const Coord end = box.max();
  int count = 0;
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        Coord c(i,j,k);
        if (!idx.isOn(c)) {
          continue;
        }
        count++;
      }  // for k
    }  // for j
  }  // for i
  archive(count);
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        Coord c(i,j,k);
        if (!idx.isOn(c)) {
          continue;
        }
        VIdx id = idx.get(c);
        assert(id < data_->size());
        archive(i,j,k,data_->at(id));
      }  // for k
    }  // for j
  }  // for i
}  // SerializeDense

template<typename T> template<class ScalarGridVIdx>
void Vector<T>::DeserializeGhost(
  int ghost_width, const ScalarGridVIdx &idx,
  canary::CanaryInputArchive &archive) {
  DeserializeHelper(idx, archive);
}  // DeserializeGhost

template<typename T> template<class ScalarGridVIdx>
void Vector<T>::DeserializeHelper(
  const ScalarGridVIdx &idx, canary::CanaryInputArchive &archive) {
  typedef typename ScalarGridVIdx::ValueType VIdx;
  int count = 0;
  archive(count);
  for (int c = 0; c < count; ++c) {
    int i=0, j=0, k=0;
    T value;
    archive(i,j,k,value);
    VIdx id = idx.get(i,j,k);
    assert(id < data_->size());
    data_->at(id) = value;
  }  // for c
}  // DeserializeHelper

namespace {
typedef int32_t VIndex;
typedef ScalarGrid<3, 3, VIndex> ScalarGridVIdx;
}  // namespace anonymous

template void Vector<float>::SendGhostData<ScalarGridVIdx>(
  canary::CanaryTaskContext*, int, const ScalarGridVIdx&) const;
template int Vector<float>::ReceiveGhostData<ScalarGridVIdx>(
  canary::CanaryTaskContext*, int, const ScalarGridVIdx&);
template void Vector<double>::SendGhostData<ScalarGridVIdx>(
  canary::CanaryTaskContext*, int, const ScalarGridVIdx&) const;
template int Vector<double>::ReceiveGhostData<ScalarGridVIdx>(
  canary::CanaryTaskContext*, int, const ScalarGridVIdx&);

}  // namespace common
