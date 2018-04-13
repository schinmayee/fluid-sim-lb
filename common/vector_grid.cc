#include <iostream>
#include <string>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include "common/vector_grid.h"

namespace common {

template<Index N1, Index N2, typename T>
VectorGrid<N1, N2, T>::VectorGrid()
	: metagrid_(), staggered_(true),
    background_(T(0)), voxel_len_(1), name_("data"),
		data_(GridT::create(background_)), accessor_{data_->getAccessor()} {
		data_->setName("data");
		data_->setGridClass(openvdb::GRID_STAGGERED);
}

// bbox is from 0,0,0 to x,y,z to include face values for staggered grid.
// When accessing velocities, should ignore final y,z for x component,
// final x,z for y component, and final x,y for z component.
template<Index N1, Index N2, typename T>
VectorGrid<N1, N2, T>::VectorGrid(
		int x, int y, int z, bool staggered, float voxel_len, Vec3T background,
		std::string name)
	: metagrid_(x,y,z), staggered_(staggered),
		background_(background), voxel_len_(voxel_len), name_(name),
		data_(GridT::create(background)), accessor_{data_->getAccessor()} {
		data_->setName(name);
		if (staggered) {
			data_->setGridClass(openvdb::GRID_STAGGERED);
		}
}

template<Index N1, Index N2, typename T>
VectorGrid<N1, N2, T>::VectorGrid(
		Coord global_dims, bool staggered, float voxel_len, Vec3T background,
		std::string name)
	: metagrid_(global_dims), staggered_(staggered),
		background_(background), voxel_len_(voxel_len), name_(name),
		data_(GridT::create(background)), accessor_{data_->getAccessor()} {
		data_->setName(name);
		if (staggered) {
			data_->setGridClass(openvdb::GRID_STAGGERED);
		}
}

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::InitializePartition(
  const MetaGrid &metagrid,
  bool staggered,
  float voxel_len, Vec3T background, std::string name) {
  metagrid_ = metagrid;
  background_ = background;
  staggered_ = staggered;
  voxel_len_ = voxel_len;
  name_ = name;
  data_ = GridT::create(background);
  accessor_ = data_->getAccessor();
  data_->setName(name);
  if (staggered) {
    data_->setGridClass(openvdb::GRID_STAGGERED);
  }
}

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::Clear() {
	data_->clear();
}

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::Fill(Vec3T value) {
  Coord start = metagrid_.start();
  Coord end = metagrid_.end() + Coord(1);
  CoordBBox bbox(start, end);
	data_->fill(bbox, value);
}

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::Prune(T threshold) {
	data_->pruneGrid(threshold);
}

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::CopyTo(VectorGrid<N1, N2, T> &to) {
	to.data_ = data_->deepCopy();
}

template<Index N1, Index N2, typename T>
Vec3<T> VectorGrid<N1, N2, T>::Interpolate(
		const Vec3<T> &ijk, int spatial_order) const {
  assert(spatial_order == 1);
	Vec3<T> result = openvdb::tools::StaggeredBoxSampler::sample(accessor_, ijk);
	return result;
}  // Interpolate

template<Index N1, Index N2, typename T>
T VectorGrid<N1, N2, T>::InterpolateComponent(
		const Vec3<T> &ijk, int component, int spatial_order) const {
  assert(spatial_order == 1);
  Vec3<T> offset(0);
  offset[component] = 0.5;
  Vec3<T> result = openvdb::tools::BoxSampler::sample(accessor_, ijk+offset);
  return result[component];
}  // InterpolateComponent

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::Load(std::string file_name) {
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
}  // Load

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::WriteVDB(std::string file_name) const {
	openvdb::GridPtrVec grids;
	openvdb::io::File file(file_name);
	grids.push_back(data_);
	file.write(grids);
	file.close();
}  // WriteVDB

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::Write(std::string file_name, bool text) const {
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
  end = end + Coord(1);  // add 1 for staggered grids
  for (typename GridT::ValueOnCIter iter = data_->cbeginValueOn();
       iter.test(); ++iter) {
    Coord c = iter.getCoord();
    Index i = c[0];
    Index j = c[1];
    Index k = c[2];
    if (i < start[0] || j < start[1] || k < start[2] ||
        i > end[0] || j > end[1] || k > end[2]) {
      continue;
    }
    Vec3T value = iter.getValue();
    fprintf(fp, "%d %d %d : %f %f %f\n",
          i, j, k, value[0], value[1], value[2]);
  }
	fclose(fp);
}  // Write

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::SendGhostData(
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
int VectorGrid<N1, N2, T>::ReceiveGhostData(
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
int VectorGrid<N1, N2, T>::SerializeGhost(
  int di, int dj, int dk, int ghost_width,
  canary::CanaryOutputArchive &archive) const {
  return SerializeGhost(Coord(di,dj,dk), ghost_width, archive);
}  // SerializeGhost

template<Index N1, Index N2, typename T>
int VectorGrid<N1, N2, T>::SerializeGhost(
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
void VectorGrid<N1, N2, T>::SerializeDense(
  const CoordBBox box, canary::CanaryOutputArchive &archive) const {
  VLOG(1) << "Serializing vector grid " << metagrid_.rank();
  archive(kMagicNumber);
  const Coord start = box.min();
  const Coord end = box.max() + Coord(1,1,1);
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
          Vec3T value = get(c);
          archive(i,j,k,value[0],value[1],value[2]);
          total++;
        }  // for k
      }  // for j
    }  // for i
    CHECK_EQ(total, count);
  }
  archive(kMagicNumber);
  VLOG(1) << "Serialized vector grid " << metagrid_.rank();
}  // SerializeDense

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::DeserializeGhost(
  int ghost_width, canary::CanaryInputArchive &archive) {
  DeserializeHelper(archive);
}  // DeserializeGhost

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::DeserializeHelper(
  canary::CanaryInputArchive &archive) {
  VLOG(1) << "Deserializing vector grid " << metagrid_.rank();
  int magic_number = 0;
  archive(magic_number);
  CHECK_EQ(magic_number, kMagicNumber);
  int count = 0;
  archive(count);
  for (int c = 0; c < count; ++c) {
    int i=0, j=0, k=0;
    T v1=0, v2=0, v3=0;
    archive(i,j,k,v1,v2,v3);
    set(i,j,k, Vec3T(v1,v2,v3));
  }  // for c
  magic_number = 0;
  archive(magic_number);
  CHECK_EQ(magic_number, kMagicNumber);
  VLOG(1) << "Deserialized vector grid " << metagrid_.rank();
}  // DeserializeHelper

template<Index N1, Index N2, typename T>
void VectorGrid<N1, N2, T>::ExtrapolateToBoundary() {
  const common::Coord start = metagrid_.start();
  const common::Coord end = metagrid_.end();
  const common::Coord global_dims = metagrid_.global_dims();
  // for loops are over end + 1 since this is a staggered mac grid =>
  // one extra value.
  // x component, lower
  if (start.x() == 0) {
    for (int j = start[1]; j <= end[1]+1; ++j) {
    for (int k = start[2]; k <= end[2]+1; ++k) {
      Vec3T v_out = get(-1,j,k);
      Vec3T v_in = get(0,j,k);
      v_out[1] = v_in[1];
      v_out[2] = v_in[2];
      set(-1,j,k,v_out);
    }  // for k
    }  // for j
  }  // if start.x() == 0
  // x component, upper
  if (end.x()+1 == global_dims[0]) {
    for (int j = start[1]; j <= end[1]+1; ++j) {
    for (int k = start[2]; k <= end[2]+1; ++k) {
      Vec3T v_out = get(global_dims[0],j,k);
      Vec3T v_in = get(global_dims[0]-1,j,k);
      v_out[1] = v_in[1];
      v_out[2] = v_in[2];
      set(global_dims[0],j,k,v_out);
    }  // for k
    }  // for j
  }  // if end.x()+1 == global_dims[0]
  // y component, lower
  if (start.y() == 0) {
    for (int i = start[0]; i <= end[0]+1; ++i) {
    for (int k = start[2]; k <= end[2]+1; ++k) {
      Vec3T v_out = get(i,-1,k);
      Vec3T v_in = get(i,0,k);
      v_out[0] = v_in[0];
      v_out[2] = v_in[2];
      set(i,-1,k,v_out);
    }  // for k
    }  // for i
  }  // if start.y() == 0
  // y component, upper
  if (end.y()+1 == global_dims[1]) {
    for (int i = start[0]; i <= end[0]+1; ++i) {
    for (int k = start[2]; k <= end[2]+1; ++k) {
      Vec3T v_out = get(i,global_dims[1],k);
      Vec3T v_in = get(i,global_dims[1]-1,k);
      v_out[0] = v_in[0];
      v_out[2] = v_in[2];
      set(i,global_dims[1],k,v_out);
    }  // for k
    }  // for i
  }  // if end.y()+1 == global_dims[1]
  // z component, lower
  if (start.z() == 0) {
    for (int i = start[0]; i <= end[0]+1; ++i) {
    for (int j = start[1]; j <= end[1]+1; ++j) {
      Vec3T v_out = get(i,j,-1);
      Vec3T v_in = get(i,j,0);
      v_out[0] = v_in[0];
      v_out[1] = v_in[1];
      set(i,j,-1,v_out);
    }  // for k
    }  // for i
  }  // if start.z() == 0
  // z component, upper
  if (end.z()+1 == global_dims[2]) {
    for (int i = start[0]; i <= end[0]+1; ++i) {
    for (int j = start[1]; j <= end[1]+1; ++j) {
      Vec3T v_out = get(i,j,global_dims[2]);
      Vec3T v_in = get(i,j,global_dims[2]-1);
      v_out[0] = v_in[0];
      v_out[1] = v_in[1];
      set(i,j,global_dims[2],v_out);
    }  // for k
    }  // for i
  }  // if end.z()+1 == 0
}  // ExtrapolateToBoundary

template class VectorGrid<3,3,float>;
template class VectorGrid<3,3,double>;

}  // namespace common
