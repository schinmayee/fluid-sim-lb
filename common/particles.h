#ifndef COMMON_PARTICLES_H
#define COMMON_PARTICLES_H

#include <cereal/archives/xml.hpp>
#include <cmath>
#include <functional>
#include <iostream>
#include <list>
#include <vector>

#include "canary/canary.h"
#include "common/metagrid.h"
#include "common/primitives.h"

namespace common {

// TODO: add physical size to measure how much defragment helps
// and to give the ability to tune frequency of calls to defragment.

// Particles class contains all particles over a partition, and includes
// metagrid information.
template<Index N1, Index N2, typename T = float>
class Particles {
  // For now, this is a vector pf particle data.
  // May want to make this a linked list of buckets per cell in a sparse grid.
	public:
		Particles();
		Particles(int x, int y, int z, T voxel_len=1, std::string name="data");
		Particles(Coord global_dims, T voxel_len=1, std::string name="data");
    ~Particles();

    // ParticleData contains position and velocity data for each particle.
    typedef struct ParticleData {
      Vec3<T> position;
      Vec3<T> velocity;
      Vec3<T> grid_velocity;
    } ParticleData;  // struct ParticleData

    // Each voxel in Particles (grid) contains the following linked list of
    // buckets of particles.
    class ParticleList {
      private:
        // Maximum number of particles to keep per bucket.
        static const int particles_per_bucket_ = 8;
        // Outer container.
        Particles *grid_;
        // Total number of particles in this list.
        int size_;  // without removed particles
        int total_size_;  // with removed particles

      public:
        ParticleList() = delete;
        ParticleList(Particles *grid) : grid_(grid), size_(0),
                                        total_size_(0) {};
        ~ParticleList() {
          if (grid_) {
            assert(size_ <= grid_->size_);
            grid_->size_ -= size_;
            assert(total_size_ <= grid_->total_size_);
            grid_->total_size_ -= total_size_;
          }
        }

        void set_grid(Particles *grid) {
          grid_ = grid;
        }  // set_grid

        // ParticleBucket has a capacity of particles_per_bucket number of particles,
        // and also stores the last used index (number of active particles).
        class ParticleBucket {
          private:
            friend class ParticleList;
            ParticleData data_[particles_per_bucket_];
            bool removed_[particles_per_bucket_];
            int next_;

          public:
            ParticleBucket() {
              for (int p = 0; p < particles_per_bucket_; ++p) {
                removed_[p] = false;
              }  // for p
              next_ = 0;
            }  // ParticleBucket
            // Returns total number of particles in the bucket, including
            // removed particles.
            inline int size() const {
              return next_;
            }  // size
            // Returns a reference to the particle.
            inline ParticleData& particle(int id) {
              if (id >= particles_per_bucket_) {
                std::cerr << "Particle id " << id << " exceeds bucket capacity"
                          << particles_per_bucket_ << std::endl;
              }
              return data_[id];
            }  // particle
            inline const ParticleData& particle(int id) const {
              if (id >= particles_per_bucket_) {
                std::cerr << "Particle id " << id << " exceeds bucket capacity"
                          << particles_per_bucket_ << std::endl;
              }
              return data_[id];
            }  // particle
            inline const ParticleData& const_particle(int id) const {
              if (id >= particles_per_bucket_) {
                std::cerr << "Particle id " << id << " exceeds bucket capacity"
                          << particles_per_bucket_ << std::endl;
              }
              return data_[id];
            }  // particle
            // Check if a particle is marked as removed.
            inline bool is_removed(int id) const {
              if (id >= particles_per_bucket_) {
                std::cerr << "Particle id " << id << " exceeds bucket capacity"
                          << particles_per_bucket_ << std::endl;
                return false;
              }
              return removed_[id];
            }  // is_removed

          private:
            inline bool AddParticle(const ParticleData &p) {
              assert(next_ < particles_per_bucket_);
              if (next_ >= particles_per_bucket_) {
                std::cerr << "Cannot add more particles to bucket" << std::endl;
                return false;
              }
              data_[next_] = p;
              removed_[next_] = false;
              next_++;
              return true;
            }  // AddParticle
            void Defragment();
        };  // class ParticleBucket

        // List of buckets containing particles,
        typedef std::list< ParticleBucket > ParticleBucketList;
        // List iterators.
        typedef typename ParticleBucketList::iterator ParticleListIterator;
        typedef typename ParticleBucketList::const_iterator ParticleListCIterator;

      private:
        // ParticleBucketList data.
        ParticleBucketList data_;

      public:
        // ParticleList iterators.
        ParticleListIterator begin() {
          return data_.begin();
        }  // begin
        ParticleListIterator end() {
          return data_.end();
        }  // end
        ParticleListCIterator cbegin() const {
          return data_.cbegin();
        }  // cbegin
        ParticleListCIterator cend() const {
          return data_.cend();
        }  // cend

        // Execute op for each non-removed particle.
        int ForEachConst(std::function<void(const ParticleData&, const Coord, int)> op,
                         const Coord c, int total=0) const {
          int count = 0;
          for (ParticleListCIterator bucket = cbegin();
               bucket != cend(); ++bucket) {
            for (int i = 0; i < bucket->size(); ++i) {
              if (bucket->is_removed(i)) {
                continue;
              }
              const ParticleData &p = bucket->const_particle(i);
              op(p, c, total);
              total++;
              count++;
            }  // go over each particle in each bucket
          }  // go over each bucket
          assert(count == size_);
          return total;
        }  // ForEachConst

        // Execute op for each non-removed particle.
        int ForEach(std::function<void(ParticleData&, const Coord, int)> op,
                    const Coord c, int total=0) {
          for (ParticleListIterator bucket = begin();
               bucket != end(); ++bucket) {
            for (int i = 0; i < bucket->size(); ++i) {
              if (bucket->is_removed(i)) {
                continue;
              }
              ParticleData &p = bucket->particle(i);
              op(p, c, total);
              total++;
            }  // go over each particle in each bucket
          }  // go over each bucket
          return total;
        }  // ForEach

        // Total number of particles in this list.
        inline int size() const { return size_; }
        inline int total_size() const { return total_size_; }

        // Marks a particle as removed, given bucket (iter) and id within the
        // bucket. Updates the size (total number of particles in the list).
        inline void MarkRemoved(ParticleListIterator iter, int id) {
          if (id >= particles_per_bucket_) {
            std::cerr << "Particle id " << id << " exceeds bucket capacity"
                      << particles_per_bucket_ << std::endl;
            return;
          }
          // If particle is not removed, remove it.
          if (iter->removed_[id] == false) {
            iter->removed_[id] = true;
            assert(size_ > 0);
            size_--;
            if (grid_) {
              grid_->size_--;
            }
          }
        }  // MarkRemoved

        void AddParticle(const ParticleData &p) {
          if (data_.empty()) {
            data_.emplace_back();
          }
          if (data_.back().next_ == particles_per_bucket_) {
            // Bucket is full, allocate new bucket.
            data_.emplace_back();
          }
          bool success = data_.back().AddParticle(p);
          assert(success);
          size_++;
          total_size_++;
          if (grid_) {
            grid_->size_++;
            grid_->total_size_++;
          }
        }  // AddParticle

        void AddParticles(const ParticleList &particles) {
          for (auto bucket = particles.cbegin(); bucket != particles.cend();
               ++bucket) {
            // Go over each particle in each bucket and add all particles that
            // are not removed.
            for (int i = 0; i < bucket->size(); ++i) {
              if (!bucket->is_removed(i)) {
                AddParticle(bucket->particle(i));
              }  // if !is_removed
            }  // for i
          }  // for bucket
        }  // AddParticles
    };  // class ParticleList

    typedef uint64_t Address;

    // Tree to store particles sparsely -- each voxel contains memory address
    // of allocated list. Voxel = 0 means no list is allocated.
		typedef typename openvdb::tree::Tree3<Address, N1, N2>::Type TreeT;
		typedef typename openvdb::Grid<TreeT> GridT;

	private:
		// Members
    // Stores partition information
    MetaGrid metagrid_;
    // Length of voxel
		T voxel_len_;
    // Data name.
		std::string name_;

    // OpenVDB grid stores pointer to list of particles.
		typename GridT::Ptr data_;
    // Buffer to store particles that have moved out to scatter particles.
    // External code can access data_, get lists at different coordinates for
    // data_, and perform operations such as add particles, mark particle as
    // removed on lists directly. But adding particles to/clearing particles
    // from buffer_ must be done with AddParticleToBuffer and ClearBuffer --
    // there is no external access to lists at different coordinates in
    // buffer_. Merging particles from buffer_ into data_ must be done with
    // MergeBufferIntoData.
    typename GridT::Ptr buffer_;

    // Accessors
		typename GridT::Accessor data_accessor_;
		typename GridT::Accessor buffer_accessor_;
    // Bounding box of active values (particles).
    CoordBBox active_bbox_;

    // This stores the total number of non-removed particles.
    int size_;
    // Size including removed particles.
    int total_size_;

	public:
    // Get metagrid.
    const MetaGrid& metagrid() const {
      return metagrid_;
    }  // metagrid()

    // Initialize partition information.
    void InitializePartition(
      const MetaGrid &metagrid, float voxel_len = 1, std::string name = "data");

    // Size, number of non-removed particles in the entire container.
    inline int size() const { return size_; }
    inline int total_size() const { return total_size_; }
    typename GridT::Ptr data() { return data_; }

    // Get particle list -- can do this only for data_, useful for iterating
    // over all particles for stepping particles for instance.
    // Don't use this to add particles. Use AddParticle/AddParticles to add
    // particles. Use this for only iterating over all particles.
    inline ParticleList* particle_list(int i, int j, int k) {
      return particle_list(Coord(i,j,k));
    }  // particle_list
    inline ParticleList* particle_list(Coord ijk) {
      if (!data_accessor_.isValueOn(ijk)) {
        return nullptr;
      }
      Address addr = data_accessor_.getValue(ijk);
      if (addr == 0) {
        return nullptr;
      }
      ParticleList* particles = reinterpret_cast<ParticleList*>(addr);
      return particles;
    }  // particle_list
    const ParticleList* const_particle_list(int i, int j, int k) const {
      return const_particle_list(Coord(i,j,k));
    }  // const_particle_list
    const ParticleList* const_particle_list(Coord ijk) const {
      Address addr = data_accessor_.getValue(ijk);
      if (addr == 0) {
        return nullptr;
      }
      const ParticleList *particles =
        reinterpret_cast<const ParticleList*>(addr);
      return particles;
    }  // const_particle_list

    // Clear all particles.
    inline void ClearData() {
      ClearParticles(data_);
    }  // ClearData
    void ClearBuffer() {
      ClearParticles(buffer_);
    }  // ClearBuffer

    // Add particle.
    inline void AddParticleToData(int i, int j, int k, const ParticleData &p) {
      AddParticle(this, data_accessor_, Coord(i,j,k), p);
    }  // AddParticleToData
    inline void AddParticleToData(Coord ijk, const ParticleData &p) {
      AddParticle(this, data_accessor_, ijk, p);
    }  // AddParticleToData
    // Add particle to buffer_, this should not affect size of particles in
    // container -- pass nullptr for outer container.
    inline void AddParticleToBuffer(int i, int j, int k, const ParticleData &p) {
      AddParticle(nullptr, buffer_accessor_, Coord(i,j,k), p);
    }  // AddParticleToBuffer
    inline void AddParticleToBuffer(Coord ijk, const ParticleData &p) {
      AddParticle(nullptr, buffer_accessor_, ijk, p);
    }  // AddParticleToBuffer

    // Merge buffer_ into data_.
    void MergeBufferIntoData();
    // Defragment data_.
    void Defragment();

    // Execute op for each particle.
    void ForEachConst(std::function<void(const ParticleData&, const Coord, int)> op) const;
    void ForEach(std::function<void(ParticleData&, const Coord, int)> op);

		// IO
		void Load(std::string file_name);
		void Write(std::string file_name, bool text=false) const;
    void WriteCheckpoint() const {
      std::string checkpoint = "/tmp/" + name_ + "_" +
                               std::to_string(metagrid_.rank());
      Write(checkpoint, false);
    }

		// Integrator steps particles using velocity given by v.
		template<typename VectorGridT, bool staggered=true,
			      int temporal=1, int spatial=1>
		void StepInGrid(VectorGridT &v, T dt);

    // Move particles to the correct cells.
    void UpdateParticleCells();

    // Update particle boundinb box.
    void EvaluateActiveBoundingBox() {
      active_bbox_ = data_->evalActiveVoxelBoundingBox();
    }

    // Delete outside particles,
    void DeleteOutsideParticles();

    // Send and receive data.
    void SendGhostData(
      canary::CanaryTaskContext *task_context, int ghost_width,
      bool send_outside_particles) const;
    int ReceiveGhostData(
      canary::CanaryTaskContext *task_context, int ghost_width);

    // Serialization and deserialization methods.
    // Serialize ghost regions and return rank for the neighboring partition.
    // If there is no neighbor in the given direction (boundary), return -1.
    int SerializeGhost(
      int di, int dj, int dk, int ghost_width,
      bool send_outside_particles, canary::CanaryOutputArchive &archive) const;
    int SerializeGhost(
      const Coord neighbor_dir, int ghost_width,
      bool send_outside_particles, canary::CanaryOutputArchive &archive) const;
    void DeserializeGhost(
      int ghost_width, canary::CanaryInputArchive &archive);

    // Method needed for archive to serialize/deserialize data when migrating.
    template <class Archive>
    void save(Archive &ar) const {
      LOG(INFO) << "Saving particles ...";
      ar(metagrid_, voxel_len_, name_);
      CoordBBox bbox = data_->evalActiveVoxelBoundingBox();
      SerializeDense<true>(bbox, ar);
      //WriteCheckpoint();
      LOG(INFO) << "Saved " << size() << " particles";
    }
    template <class Archive>
    void load(Archive &ar) {
      LOG(INFO) << "Loading particles ...";
      ar(metagrid_, voxel_len_, name_);
      InitializePartition(metagrid_, voxel_len_, name_);
      DeserializeHelper<true>(ar);
      EvaluateActiveBoundingBox();
      LOG(INFO) << "Loaded " << size() << " particles";
    }

  private:
    // Helper method -- free all allocated memory.
    static void ClearParticles(typename GridT::Ptr data);

    // Helper method -- add a particle.
    inline static void AddParticle(
      Particles *grid, typename GridT::Accessor &data_accessor,
      Coord ijk, const ParticleData &p) {
      Address addr = data_accessor.getValue(ijk);
      ParticleList *particles = nullptr;
      if (addr == 0) {
        particles = new ParticleList(grid);
        data_accessor.setValue(ijk, reinterpret_cast<Address>(particles));
      } else {
        particles = reinterpret_cast<ParticleList*>(addr);
      }
      particles->AddParticle(p);
    }  // AddParticle

    // Helper method -- add all particles from a list.
    inline static void AddParticles(
      Particles *grid, typename GridT::Accessor &data_accessor,
      Coord ijk, const ParticleList &add_list) {
      Address addr = data_accessor.getValue(ijk);
      ParticleList *particles = nullptr;
      if (addr == 0) {
        particles = new ParticleList(grid);
        data_accessor.setValue(ijk, reinterpret_cast<Address>(particles));
      } else {
        particles = reinterpret_cast<ParticleList*>(addr);
      }
      particles->AddParticles(add_list);
    }  // AddParticles

    inline static void WriteParticle(FILE *fp, const ParticleData &p,
                                     const Coord c, int id) {
      fprintf(fp, "%d : %f %f %f, %f %f %f\n", id,
          p.position[0], p.position[1], p.position[2],
          p.velocity[0], p.velocity[1], p.velocity[2]);
    }  // WriteParticle

    // Helper serialization/deserialization methods.
    template<bool send_grid_velocity>
    void SerializeDense(
      const CoordBBox box, canary::CanaryOutputArchive &archive) const;
    template<bool send_grid_velocity>
    inline static void SerializeParticle(
      canary::CanaryOutputArchive &archive,
      const ParticleData &p, const Coord c, int id) {
        archive(p.position[0], p.position[1], p.position[2],
                p.velocity[0], p.velocity[1], p.velocity[2]);
        if (send_grid_velocity) {
          archive(p.grid_velocity[0], p.grid_velocity[1], p.grid_velocity[2]);
        }
    }
    template<bool receive_grid_velocity>
    void DeserializeHelper(canary::CanaryInputArchive &archive);

    // Count particles in a region, for serializing.
    int Count(const CoordBBox box) const;

		//// Integrator steps particles using velocity given by v.
		//template<typename VelocityIntegrator, int temporal>
		//inline static void StepParticleInGrid(
    //  // Total number of particles in this list.
    //  VelocityIntegrator &v, T dt, ParticleData &p, const Coord c, int id) {
    //  p.position += p.grid_velocity*dt;
    //  //typedef typename VelocityIntegrator::VecType VecT;
    //  //v.template rungeKutta<temporal, VecT>(dt, p.position);
    //  assert(!std::isnan(p.position[0]));
    //  assert(!std::isnan(p.position[1]));
    //  assert(!std::isnan(p.position[2]));
    //}  // StepParticleInGrid

};  // class Particles

}  // namespace common

#endif  // COMMON_PARTICLES_H
