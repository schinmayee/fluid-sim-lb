#include <iostream>
#include <string>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include "common/scalar_grid.h"

namespace common {

template<Index N1, Index N2, typename T>
ScalarGrid<N1, N2, T>::ScalarGrid()
  : metagrid_(), background_(T(0)), voxel_len_(1), name_("data"),
    data_(GridT::create(background_)), accessor_{data_->getAccessor()} {
  data_->setName("data");
}

template<Index N1, Index N2, typename T>
ScalarGrid<N1, N2, T>::ScalarGrid(
		int x, int y, int z, float voxel_len, T background,
		std::string name)
	: metagrid_(x,y,z), background_(background), voxel_len_(voxel_len), name_(name),
		data_(GridT::create(background)), accessor_{data_->getAccessor()} {
		data_->setName(name);
}

template<Index N1, Index N2, typename T>
ScalarGrid<N1, N2, T>::ScalarGrid(
		Coord global_dims, float voxel_len, T background, std::string name)
	: metagrid_(global_dims), background_(background), voxel_len_(voxel_len),
    name_(name),
    data_(GridT::create(background)), accessor_{data_->getAccessor()} {
		data_->setName(name);
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::InitializePartition(
  const MetaGrid &metagrid, float voxel_len, T background, std::string name) {
  metagrid_ = metagrid;
  background_ = background;
  voxel_len_ = voxel_len;
  name_ = name;
  data_ = GridT::create(background);
  accessor_ = data_->getAccessor();
  data_->setName(name);
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::Clear() {
	data_->clear();
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::Fill(T value) {
  CoordBBox bbox(metagrid_.start(), metagrid_.end());
	data_->fill(bbox, value);
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::Prune(T threshold) {
	data_->pruneGrid(threshold);
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::CopyTo(ScalarGrid<N1, N2, T> &to) {
	to.data_ = data_->deepCopy();
}

template<Index N1, Index N2, typename T>
int ScalarGrid<N1, N2, T>::Count(CoordBBox box) const {
  const Coord start = box.min();
  const Coord end = box.max();
  int total = 0;
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        if (isOn(i,j,k)) {
          total ++;
        }
      }  // for k
    }  // for j
  }  // for i
  return total;
}

template<Index N1, Index N2, typename T>
int ScalarGrid<N1, N2, T>::Count(CoordBBox box, T value) const {
  const Coord start = box.min();
  const Coord end = box.max();
  int total = 0;
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        if (get(i,j,k) == value) {
          total ++;
        }
      }  // for k
    }  // for j
  }  // for i
  return total;
}

template<Index N1, Index N2, typename T>
T ScalarGrid<N1, N2, T>::Interpolate(
		const Vec3<T> &ijk, int spatial_order) const {
  assert(spatial_order == 1);
	T result = openvdb::tools::BoxSampler::sample(accessor_, ijk);
	return result;
}  // Interpolate

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::Load(std::string file_name) {
	openvdb::io::File file(file_name);
	file.open();
	openvdb::GridBase::Ptr base_grid;
	for (openvdb::io::File::NameIterator name_iter = file.beginName();
			name_iter != file.endName(); ++name_iter) {
		if (name_iter.gridName() == name_) {
			base_grid = file.readGrid(name_iter.gridName());
		}
	}
	data_ = openvdb::gridPtrCast<GridT>(base_grid);
	accessor_ = data_->getAccessor();
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::WriteVDB(std::string file_name) const {
	openvdb::GridPtrVec grids;
	openvdb::io::File file(file_name);
	grids.push_back(data_);
	file.write(grids);
	file.close();
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::Write(std::string file_name, bool text) const {
	std::string mode = text? "wt" : "wb";
	FILE *fp = fopen(file_name.c_str(), mode.c_str());
	if (fp == NULL) {
		fprintf(stderr, "Could not create %s file\n", file_name.c_str());
		exit(-1);
	}
  Coord start = metagrid_.start();
  Coord end = metagrid_.end();
	fprintf(fp, "%d %d %d\n", start[0], start[1], start[2]);
	fprintf(fp, "%d %d %d\n", end[0], end[1], end[2]);
  for (typename GridT::ValueOnCIter iter = data_->cbeginValueOn(); iter.test(); ++iter) {
    Coord c = iter.getCoord();
    Index i = c[0];
    Index j = c[1];
    Index k = c[2];
    if (i < start[0] || j < start[1] || k < start[2] ||
        i > end[0] || j > end[1] || k > end[2]) {
      continue;
    }
    float value = iter.getValue();
		fprintf(fp, "%d %d %d : %f\n", i, j, k, value);
  }
	fclose(fp);
}

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::SendGhostData(
  canary::CanaryTaskContext *task_context, int ghost_width) const {
  int count = 0;
  // Get neighbor rank and serialize data for each neighbor.
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        if (i == 0 && j == 0 && k == 0) {
          continue;
        }
        if (metagrid_.IsValidNeighbor(i,j,k)) {
          struct evbuffer *buffer = evbuffer_new();
          canary::CanaryOutputArchive archive(buffer);
          int rank = SerializeGhost(
            i, j, k, ghost_width, archive);
          assert(rank >= 0);
          task_context->Scatter(rank, canary::RawEvbuffer(buffer));
          count++;
        }  // if valid neighbor
      }  // for k
    }  // for j
  }  // for i
  assert(count == metagrid_.NeighborCount());
}  // SendGhostData

template<Index N1, Index N2, typename T>
int ScalarGrid<N1, N2, T>::ReceiveGhostData(
  canary::CanaryTaskContext *task_context, int ghost_width) {
  // Go over data rceived in receive buffer.
  EXPECT_GATHER_SIZE(metagrid_.NeighborCount());
  auto recv_buffer = task_context->Gather<canary::RawEvbuffer>();
  for (auto &raw_buffer : recv_buffer) {
    canary::CanaryInputArchive archive(raw_buffer.buffer);
    DeserializeGhost(ghost_width, archive);
    evbuffer_free(raw_buffer.buffer);
  }  // for raw_buffer
  return 0;
}  // ReceiveGhostData

template<Index N1, Index N2, typename T>
int ScalarGrid<N1, N2, T>::SerializeGhost(
  int di, int dj, int dk, int ghost_width,
  canary::CanaryOutputArchive &archive) const {
  return SerializeGhost(Coord(di,dj,dk), ghost_width, archive);
}  // SerializeGhost

template<Index N1, Index N2, typename T>
int ScalarGrid<N1, N2, T>::SerializeGhost(
  Coord neighbor_dir, int ghost_width,
  canary::CanaryOutputArchive &archive) const {
  const Coord pid = metagrid_.partition_id();
  const Coord neighbor_id = pid + neighbor_dir;
  if (!metagrid_.IsValidPartitionID(neighbor_id)) {
    VLOG(1) << "No neighbor in " << ToString(neighbor_dir) << " for "
              << ToString(pid);
    return -1;
  }
  CoordBBox bbox = metagrid_.ComputeBBox(neighbor_id);
  VLOG(1) << "Original box for neighbor " << ToString(neighbor_dir) <<
               " of " << ToString(metagrid_.partition_id()) << " is " <<
               ToString(bbox);
  bbox.expand(ghost_width);
  bbox.intersect(metagrid_.bbox());
  VLOG(1) << "Box for neighbor " << ToString(neighbor_dir) << " of " <<
               ToString(metagrid_.partition_id()) << " is " <<
               ToString(bbox); 
  SerializeDense(bbox, archive);
  int rank = metagrid_.ComputeRank(neighbor_id);
  return rank;
}  // SerializeGhost

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::SerializeDense(
  const CoordBBox box, canary::CanaryOutputArchive &archive) const {
  VLOG(1) << "Serializing scalar grid " << metagrid_.rank();
  archive(kMagicNumber);
  const Coord start = box.min();
  const Coord end = box.max();
  int count = 0;
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        Coord c(i,j,k);
        if (!isOn(c)) {
          continue;
        }
        count++;
      }  // for k
    }  // for j
  }  // for i
  archive(count);
  if (count != 0) {
    int total = 0;
    for (int i = start[0]; i <= end[0]; ++i) {
      for (int j = start[1]; j <= end[1]; ++j) {
        for (int k = start[2]; k <= end[2]; ++k) {
          Coord c(i,j,k);
          if (!isOn(c)) {
            continue;
          }
          archive(i,j,k,get(c));
          total++;
        }  // for k
      }  // for j
    }  // for i
    CHECK_EQ(total, count);
  }
  archive(kMagicNumber);
  VLOG(1) << "Serialized scalar grid " << metagrid_.rank();
}  // SerializeDense

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::DeserializeGhost(
  int ghost_width, canary::CanaryInputArchive &archive) {
  DeserializeHelper(archive);
}  // DeserializeGhost

template<Index N1, Index N2, typename T>
void ScalarGrid<N1, N2, T>::DeserializeHelper(
  canary::CanaryInputArchive &archive) {
  VLOG(1) << "Deserializing scalar grid " << metagrid_.rank();
  int magic_number = 0;
  archive(magic_number);
  CHECK_EQ(magic_number, kMagicNumber);
  int count = 0;
  archive(count);
  for (int c = 0; c < count; ++c) {
    int i=0, j=0, k=0;
    T value;
    archive(i,j,k,value);
    set(i,j,k, value);
  }  // for c
  magic_number = 0;
  archive(magic_number);
  CHECK_EQ(magic_number, kMagicNumber);
  VLOG(1) << "Deserialized scalar grid " << metagrid_.rank();
}  // DeserializeHelper

template class ScalarGrid<3,3,float>;
template class ScalarGrid<3,3,double>;
template class ScalarGrid<3,3,int64_t>;
template class ScalarGrid<3,3,int32_t>;
template class ScalarGrid<3,3,bool>;

}  // namespace common
