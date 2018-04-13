#include <assert.h>
#include <cmath>
#include <iostream>
#include <openvdb/tools/VelocityFields.h>

#include "common/particles.h"
#include "common/vector_grid.h"

namespace common {

template<Index N1, Index N2, typename T>
Particles<N1, N2, T>::Particles()
  : metagrid_(), voxel_len_(1), name_("data"),
    data_(GridT::create(Address(0))), buffer_(GridT::create(Address(0))),
    data_accessor_(data_->getAccessor()),
    buffer_accessor_(buffer_->getAccessor()),
    size_(0), total_size_(0) {
  data_->setName("data");
}  // Particles

template<Index N1, Index N2, typename T>
Particles<N1, N2, T>::Particles(
		int x, int y, int z, T voxel_len,
		std::string name)
  : metagrid_(x,y,z), voxel_len_(voxel_len), name_(name),
    data_(GridT::create(Address(0))), buffer_(GridT::create(Address(0))),
    data_accessor_(data_->getAccessor()),
    buffer_accessor_(buffer_->getAccessor()),
    size_(0), total_size_(0) {
  data_->setName(name);
}  // Particles

template<Index N1, Index N2, typename T>
Particles<N1, N2, T>::Particles(
		Coord global_dims, T voxel_len,
		std::string name)
  : metagrid_(global_dims), voxel_len_(voxel_len), name_(name),
    data_(GridT::create(Address(0))), buffer_(GridT::create(Address(0))),
    data_accessor_(data_->getAccessor()),
    buffer_accessor_(buffer_->getAccessor()),
    size_(0), total_size_(0) {
  data_->setName("data");
}  // Particles

template<Index N1, Index N2, typename T>
Particles<N1, N2, T>::~Particles() {
  // Free all allocated memory.
  ClearData();
  ClearBuffer();
  size_ = 0;
  total_size_ = 0;
}  // ~Particles

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::InitializePartition(
  const MetaGrid &metagrid, float voxel_len, std::string name) {
  metagrid_ = metagrid;
  voxel_len_ = voxel_len;
  name_ = name;
  ClearData();
  ClearBuffer();
  size_ = 0;
  total_size_ = 0;
  data_ = GridT::create(Address(0));
  buffer_ = GridT::create(Address(0));
  data_accessor_ = data_->getAccessor();
  buffer_accessor_ = buffer_->getAccessor();
  data_->setName(name);
}  // IntiializePartition

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::ClearParticles(typename GridT::Ptr data) {
  for (typename GridT::ValueOnIter iter = data->beginValueOn();
       iter.test(); ++iter) {
    Address addr = iter.getValue();
    if (addr == 0) {
      continue;
    }  // if addr == 0
    ParticleList *list = reinterpret_cast<ParticleList*>(addr);
    delete list;
    iter.setValue(0);
  }  // outermost for loop to go over active values
  data->clear();
}  // ClearParticles

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::MergeBufferIntoData() {
  for (typename GridT::ValueOnIter iter = buffer_->beginValueOn();
       iter.test(); ++iter) {
    Address addr = iter.getValue();
    if (addr == 0) {
      continue;
    }  // if addr == 0
    ParticleList *list = reinterpret_cast<ParticleList*>(addr);
    AddParticles(this, data_accessor_, iter.getCoord(), *list);
  }  // outermost for loop to go over active values
}  // MergeBufferIntoData

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::Defragment() {
  for (typename GridT::ValueOnIter iter = data_->beginValueOn();
       iter.test(); ++iter) {
    Address addr = iter.getValue();
    if (addr == 0) {
      continue;
    }  // if addr == 0
    ParticleList *list = reinterpret_cast<ParticleList*>(addr);
    if (list->size() == 0) {
      // No particles in this list.
      iter.setValue(0);
      iter.setValueOff();
    } else {
      ParticleList *defrag = new ParticleList(this);
      // Copy all particles that are not removed.
      defrag->AddParticles(*list);
      // Set new list, free old list.
      iter.setValue(reinterpret_cast<Address>(defrag));
    }
    // Free memory and update particle count in container.
    delete list;
  }  // outermost for loop to go over active voxels containing particles
  // Make the grid sparse.
  data_->pruneGrid();
}  // Defragment

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::Load(std::string file_name) {
	std::cerr << "Load particles from file unimplemented\n";
	exit(-1);
}  // Load

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::ForEachConst(
  std::function<void(const ParticleData&, const Coord, int)> op) const {
  int total = 0;
  for (typename GridT::ValueOnCIter iter = data_->cbeginValueOn();
       iter.test(); ++iter) {
    Address addr = iter.getValue();
    if (addr == 0) {
      continue;
    }  // if addr == 0
    ParticleList *list = reinterpret_cast<ParticleList*>(addr);
    total = list->ForEachConst(op, iter.getCoord(), total);
	}  // iter over active voxels containing particles
  assert(total == size_);
}  // ForEachConst

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::ForEach(
  std::function<void(ParticleData&, const Coord, int)> op) {
  int total = 0;
  for (typename GridT::ValueOnIter iter = data_->beginValueOn();
       iter.test(); ++iter) {
    Address addr = iter.getValue();
    if (addr == 0) {
      continue;
    }  // if addr == 0
    ParticleList *list = reinterpret_cast<ParticleList*>(addr);
    total = list->ForEach(op, iter.getCoord(), total);
	}  // iter over active voxels containing particles
  assert(total == size_);
}  // ForEach

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::Write(std::string file_name, bool text) const {
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
	fprintf(fp, "%d\n", size_);
  // Bind arguments and call WriteParticle for each particle.
  std::function<void(const ParticleData&, const Coord, int)> write_op =
    std::bind(WriteParticle, std::ref(fp),
              std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3);
  ForEachConst(write_op);
	fclose(fp);
}  // Write

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::DeleteOutsideParticles() {
  Coord start = metagrid_.start();
  Coord end = metagrid_.end();
  for (typename GridT::ValueOnIter iter = data_->beginValueOn();
       iter.test(); ++iter) {
    Coord c = iter.getCoord();
    if (c[0] < start[0] || c[1] < start[1] || c[2] < start[2] ||
        c[0] > end[0]   || c[1] > end[1]   || c[2] > end[2]) {
      Address addr = iter.getValue();
      if (addr == 0) {
        continue;
      }
      ParticleList *list = reinterpret_cast<ParticleList*>(addr);
      delete list;
      iter.setValue(0);
      iter.setValueOff();
    }
  }  // for iter
  data_->pruneGrid();
}  // DeleteOutsideParticles

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::SendGhostData(
  canary::CanaryTaskContext *task_context, int ghost_width,
  bool send_outside_particles) const {
  int count = 0;
  // Get neighbor rank and serialized data for each neighbor.
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <=1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        if (i == 0 && j == 0 && k == 0) {
          continue;
        }
        if (metagrid_.IsValidNeighbor(i,j,k)) {
          struct evbuffer *buffer = evbuffer_new();
          canary::CanaryOutputArchive archive(buffer);
          int rank = SerializeGhost(
            i, j, k, ghost_width, send_outside_particles, archive);
          assert(rank >= 0);
          VLOG(1) << "Scattering particles to " << rank;
          task_context->Scatter(rank, canary::RawEvbuffer(buffer));
          count++;
        }  // if valid neighbor
      }  // for k
    }  // for j
  }  // for i
  assert(count == metagrid_.NeighborCount());
}  // SendGhostData

template<Index N1, Index N2, typename T>
int Particles<N1, N2, T>::ReceiveGhostData(
  canary::CanaryTaskContext *task_context, int ghost_width) {
  // Go over data received in receive buffer.
  EXPECT_GATHER_SIZE(metagrid_.NeighborCount());
  auto recv_buffer = task_context->Gather<canary::RawEvbuffer>();
  for (auto &raw_buffer : recv_buffer) {
    canary::CanaryInputArchive archive(raw_buffer.buffer);
    VLOG(1) << "Derializing some particles ..."; 
    DeserializeGhost(ghost_width, archive);
    evbuffer_free(raw_buffer.buffer);
  }  // for raw_buffer
  return 0;
}  // ReceiveGhostData

template<Index N1, Index N2, typename T>
int Particles<N1, N2, T>::SerializeGhost(
  int di, int dj, int dk, int ghost_width,
  bool send_outside_particles, canary::CanaryOutputArchive &archive) const {
  return SerializeGhost(Coord(di,dj,dk), ghost_width,
                        send_outside_particles, archive);
}  // SerializeGhost

template<Index N1, Index N2, typename T>
int Particles<N1, N2, T>::SerializeGhost(
  const Coord neighbor_dir, int ghost_width,
  bool send_outside_particles, canary::CanaryOutputArchive &archive) const {
  Coord neighbor_id = metagrid_.partition_id() + neighbor_dir;
  if (!metagrid_.IsValidPartitionID(neighbor_id)) {
    VLOG(1) << "No neighbor in " << ToString(neighbor_dir) << " for "
              << ToString(metagrid_.partition_id());
    return -1;
  }
  CoordBBox bbox = metagrid_.ComputeBBox(neighbor_id);
  VLOG(1) << "Original box for neighbor " << ToString(neighbor_dir) <<
               " of " << ToString(metagrid_.partition_id()) << " is " <<
               ToString(bbox);
  bbox.expand(ghost_width);
  if (send_outside_particles) {
    bbox.intersect(active_bbox_);
  } else {
    bbox.intersect(active_bbox_);
    bbox.intersect(metagrid_.bbox());
  }
  VLOG(1) << "Box for neighbor " << ToString(neighbor_dir) << " of " <<
               ToString(metagrid_.partition_id()) << " is " <<
               ToString(bbox); 
  SerializeDense<false>(bbox, archive);
  int rank = metagrid_.ComputeRank(neighbor_id);
  return rank;
}  // SerializeGhost

template<Index N1, Index N2, typename T>
template<bool send_grid_velocity>
void Particles<N1, N2, T>::SerializeDense(
  const CoordBBox box, canary::CanaryOutputArchive &archive) const {
  VLOG(1) << "Serializing particles " << metagrid_.rank();
  archive(kMagicNumber);
  int count = Count(box);
  archive(count);
  const Coord start = box.min();
  const Coord end = box.max();
  if (count != 0) {
    int total = 0;
    for (int i = start[0]; i <= end[0]; ++i) {
      for (int j = start[1]; j <= end[1]; ++j) {
        for (int k = start[2]; k <= end[2]; ++k) {
          const ParticleList *list = const_particle_list(i, j, k);
          if (list == nullptr) {
            continue;
          }
          int list_size = list->size();
          if (list_size == 0) {
            continue;
          }
          std::function<void(const ParticleData&, const Coord, int)> ser_op =
            std::bind(
              SerializeParticle<send_grid_velocity>,
              std::ref(archive), std::placeholders::_1,
              std::placeholders::_2, std::placeholders::_3);
          total = list->ForEachConst(ser_op, Coord(i,j,k), total);
        }  // for k
      }  // for j
    }  // for i
    CHECK_EQ(count, total);
  }
  archive(kMagicNumber);
  // End, 0 more particles
  VLOG(1) << "Serialized particles " << metagrid_.rank();
  VLOG(1) << "Serialized " << count << " particles for " <<
               ToString(box) << " from " << ToString(metagrid_.bbox());
}  // SerializeDense

template<Index N1, Index N2, typename T>
int Particles<N1, N2, T>::Count(const CoordBBox box) const {
  const Coord start = box.min();
  const Coord end = box.max();
  int total = 0;
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        const ParticleList *list = const_particle_list(i, j, k);
        if (list == nullptr) {
          continue;
        }
        total += list->size();
      }  // for k
    }  // for j
  }  // for i
  return total;
}  // Count

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::DeserializeGhost(
  int ghost_width, canary::CanaryInputArchive &archive) {
  DeserializeHelper<false>(archive);
}  // DeserializeGhost

template<Index N1, Index N2, typename T>
template<bool receive_grid_velocity>
void Particles<N1, N2, T>::DeserializeHelper(
  canary::CanaryInputArchive &archive) {
  VLOG(1) << "Deserializing particles " << metagrid_.rank();
  int magic_number = 0;
  archive(magic_number);
  CHECK_EQ(magic_number, kMagicNumber);
  int count = 0;
  archive(count);
  for (int i = 0; i < count; ++i) {
    ParticleData p;
    archive(p.position[0], p.position[1], p.position[2],
            p.velocity[0], p.velocity[1], p.velocity[2]); 
    if (receive_grid_velocity) {
      archive(p.grid_velocity[0], p.grid_velocity[1], p.grid_velocity[2]);
    } else {
      p.grid_velocity = Vec3<T>(0);
    }
    Coord c;
    Vec3<T> pos = p.position;
    if (receive_grid_velocity) {
      pos /= voxel_len_;
    }
    c[0] = floor(pos[0]);
    c[1] = floor(pos[1]);
    c[2] = floor(pos[2]);
    AddParticleToData(c, p);
  }  // list_s
  magic_number = 0;
  archive(magic_number);
  CHECK_EQ(magic_number, kMagicNumber);
  VLOG(1) << "Deserialized particles " << metagrid_.rank();
  VLOG(1) << "Partition " << ToString(metagrid_.bbox()) << " deserialized "
            << count;
}  // DeserializeHelper

template<Index N1, Index N2, typename T>
template<typename VectorGridT, bool staggered, int temporal, int spatial>
void Particles<N1, N2, T>::StepInGrid(VectorGridT &v, T dt) {
	if (staggered) {
		assert(v.isStaggered());
	}
	// VelocityIntegrator to handle integration
	//typedef typename openvdb::tools::VelocityIntegrator<typename VectorGridT::GridT, staggered, spatial> VelocityIntegratorT;
  //VelocityIntegratorT integrator(*(v.data()));
  // Bind arguments and call StepParticleInGrid for each particle.
  //std::function<void(ParticleData&, const Coord, int)> step_op =
  //  std::bind(StepParticleInGrid<VelocityIntegratorT, temporal>,
  //            std::ref(integrator), dt,
  //            std::placeholders::_1, std::placeholders::_2,
  //            std::placeholders::_3);
  //ForEach(step_op);

  // Updateparticle cells and do defragmentation all in one step.

  // Update positions and insert particles into buffer, and clear particles.
  const Coord dims = metagrid_.global_dims();
  int num_particles = size();
  ParticleData *buffer = new ParticleData[num_particles];
  int p = 0;
  for (typename GridT::ValueOnIter iter = data_->beginValueOn();
       iter.test(); ++iter) {
    Address addr = iter.getValue();
    if (addr == 0) {
      continue;
    }  // if addr == 0
    ParticleList *list = reinterpret_cast<ParticleList*>(addr);
    for (typename ParticleList::ParticleListIterator bucket = list->begin();
         bucket != list->end(); ++bucket) {
      for (int i = 0; i < bucket->size(); ++i) {
        assert(!bucket->is_removed(i));
        //if (bucket->is_removed(i)) {
        //  continue;
        //}
        CHECK(p < num_particles);
        buffer[p] = bucket->particle(i);
        buffer[p].position += buffer[p].grid_velocity * dt;
        p++;
      }  // for i, for each particle in bucket
    }  // for bucket
    delete list;
  }  // iter over active voxels containing particles
  data_->clear();

  // Now copy particles over, copy only those particles that are still within the grid.
  for (int i = 0; i < num_particles; ++i) {
    ParticleData &pd = buffer[i];
    Coord id(0);
    for (int d = 0; d < 3; ++d) {
      //id[d] = floor(pd.position[d]/voxel_len_);
      id[d] = floor(pd.position[d]);
    }  // for d
    if (id[0] < -1 || id[0] > dims[0]+1 ||
        id[1] < -1 || id[1] > dims[1]+1 ||
        id[2] < -1 || id[2] > dims[2]+1) {
      continue;
    }
    if (id[0] < 0) {
      id[0] = 0;
      pd.position[0] = 0;
    }
    if (id[1] < 0) {
      id[1] = 0;
      pd.position[1] = 0;
    }
    if (id[2] < 0) {
      id[2] = 0;
      pd.position[2] = 0;
    }
    if (id[0] >= dims[0]) {
      id[0] = dims[0];
      pd.position[0] = T(dims[0]);
    }
    if (id[1] >= dims[1]) {
      id[1] = dims[1];
      pd.position[1] = T(dims[1]);
    }
    if (id[2] >= dims[2]) {
      id[2] = dims[2];
      pd.position[2] = T(dims[2]);
    }
    AddParticleToData(id, pd);
  }
  delete buffer;
}  // StepInGrid

template<Index N1, Index N2, typename T>
void Particles<N1, N2, T>::UpdateParticleCells() {
  ClearBuffer();
  int moved = 0;
  // VLOG(1) << "Original number of particles = " << size();
  for (typename GridT::ValueOnIter iter = data_->beginValueOn();
       iter.test(); ++iter) {
    Address addr = iter.getValue();
    if (addr == 0) {
      continue;
    }  // if addr == 0
    ParticleList *list = reinterpret_cast<ParticleList*>(addr);
    for (typename ParticleList::ParticleListIterator bucket = list->begin();
         bucket != list->end(); ++bucket) {
      for (int i = 0; i < bucket->size(); ++i) {
        if (bucket->is_removed(i)) {
          continue;
        }
        ParticleData &pd = bucket->particle(i);
        Coord id(0);
        for (int d = 0; d < 3; ++d) {
          //id[d] = floor(pd.position[d]/voxel_len_);
          id[d] = floor(pd.position[d]);
        }  // for d
        if (id != iter.getCoord()) {
          // Remove particle and add it to correct cell in buffer.
          list->MarkRemoved(bucket, i);
          AddParticleToBuffer(id, pd);
          moved++;
        }  // id != iter.getCoord()
      }  // for i, for each particle in bucket
    }  // for bucket
  }  // iter over active voxels containing particles
  // VLOG(1) << "Number of particles to be moved = " << moved;
  int old_size = size_;
  // VLOG(1) << "Number of particles before merging = " << size();
  MergeBufferIntoData();
  // VLOG(1) << "Number of particles after merging = " << size();
  assert(old_size + moved == size_);
}  // UpdateParticleCells

template class Particles<3, 3, float>;

template void Particles<3, 3, float>::SerializeDense<true>(const CoordBBox box, canary::CanaryOutputArchive&) const;
template void Particles<3, 3, float>::SerializeDense<false>(const CoordBBox box, canary::CanaryOutputArchive&) const;
template void Particles<3, 3, float>::DeserializeHelper<true>(canary::CanaryInputArchive&);
template void Particles<3, 3, float>::DeserializeHelper<false>(canary::CanaryInputArchive&);

template void Particles<3, 3, float>::StepInGrid< VectorGrid<3,3,float>, true, 1, 1 >(VectorGrid<3,3,float>&, float);
template void Particles<3, 3, float>::StepInGrid< VectorGrid<3,3,float>, true, 2, 1 >(VectorGrid<3,3,float>&, float);
template void Particles<3, 3, float>::StepInGrid< VectorGrid<3,3,float>, true, 2, 2 >(VectorGrid<3,3,float>&, float);
template void Particles<3, 3, float>::StepInGrid< VectorGrid<3,3,float>, false, 1, 1 >(VectorGrid<3,3,float>&, float);
template void Particles<3, 3, float>::StepInGrid< VectorGrid<3,3,float>, false, 2, 1 >(VectorGrid<3,3,float>&, float);
template void Particles<3, 3, float>::StepInGrid< VectorGrid<3,3,float>, false, 2, 2 >(VectorGrid<3,3,float>&, float);
template void Particles<3, 3, double>::StepInGrid< VectorGrid<3,3,double>, true, 1, 1 >(VectorGrid<3,3,double>&, double);
template void Particles<3, 3, double>::StepInGrid< VectorGrid<3,3,double>, true, 2, 1 >(VectorGrid<3,3,double>&, double);
template void Particles<3, 3, double>::StepInGrid< VectorGrid<3,3,double>, true, 2, 2 >(VectorGrid<3,3,double>&, double);
template void Particles<3, 3, double>::StepInGrid< VectorGrid<3,3,double>, false, 1, 1 >(VectorGrid<3,3,double>&, double);
template void Particles<3, 3, double>::StepInGrid< VectorGrid<3,3,double>, false, 2, 1 >(VectorGrid<3,3,double>&, double);
template void Particles<3, 3, double>::StepInGrid< VectorGrid<3,3,double>, false, 2, 2 >(VectorGrid<3,3,double>&, double);

}  // namespace common
