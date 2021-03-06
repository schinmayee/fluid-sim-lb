#ifndef COMMON_VECTOR_GRID_H
#define COMMON_VECTOR_GRID_H

#include <cereal/archives/xml.hpp>
#include <iostream>
#include <string>

#include <openvdb/openvdb.h>

#include "canary/canary.h"
#include "common/metagrid.h"
#include "common/primitives.h"

namespace common {

template<Index N1, Index N2, typename T>
class VectorGrid {
	public:
		typedef Vec3<T> Vec3T;
		typedef typename openvdb::tree::Tree3<Vec3T, N1, N2>::Type TreeT;
		typedef typename openvdb::Grid<TreeT> GridT;
    typedef T ValueType;
    template<typename ValueT>
    struct ValueConverter {
      using Type = VectorGrid<N1, N2, ValueT>;
    };

		// constructor
		VectorGrid();
		VectorGrid(int x, int y, int z,
				       bool staggered = true, float voxel_len = 1,
							 Vec3T background = Vec3T(0), std::string name = "data");
		VectorGrid(Coord global_dims, bool staggered = true, float voxel_len = 1,
							 Vec3T background = Vec3T(0), std::string name = "data");

	private:
		// members
    MetaGrid metagrid_;
		bool staggered_;  // is this a staggered grid?
		Vec3T background_;  // background value
		float voxel_len_;  // length of each voxel
		std::string name_;
		typename GridT::Ptr data_;  // openvdb grid with data
		typename GridT::Accessor accessor_;  // grid accessor

	public:
		// methods

    void InitializePartition(
      const MetaGrid &metagrid,
      bool staggered = true,
      float voxel_len = 1, Vec3T background = Vec3T(0), std::string name = "data");

		inline int dims(size_t i) const { return metagrid_.dims(i); }
		inline bool isStaggered() const { return staggered_; }
		inline float voxel_len() const { return voxel_len_; }
		inline Vec3T get(int i, int j, int k) const {
			return accessor_.getValue(Coord(i,j,k));
		}
		inline Vec3T get(Coord c) const {
			return accessor_.getValue(c);
		}
		inline void set(int i, int j, int k, Vec3T value) {
			accessor_.setValue(Coord(i,j,k), value);
		}
		inline void set(Coord c, Vec3T value) {
			accessor_.setValue(c, value);
		}
		inline bool isOn(int i, int j, int k) const {
			return accessor_.isValueOn(Coord(i,j,k));
		}
		inline bool isOn(Coord c) const {
			return accessor_.isValueOn(c);
		}
    inline void fill(const CoordBBox &box, Vec3T value, bool active=true) {
      accessor_.getTree()->fill(box, value, active);
    }
		inline typename GridT::Ptr data() { return data_; }
		inline typename GridT::Ptr const_data() const { return data_; }

    // Get metagrid.
    const MetaGrid& metagrid() const {
      return metagrid_;
    }  // metagrid()

		// Empty the grid, so that all voxels become inactive bakground.
		void Clear();
		// Fill all values to a given constant -- may not be optimally sparse.
		void Fill(Vec3T value);
		// Make the grid sparse by pruning.
		void Prune(T tolerance = T(0));
		// Copy to another grid.
		void CopyTo(VectorGrid<N1, N2, T> &to);

		// Interpolate to give posiition. Use spatial_order interpolation
		// (spatial_order = 1 means linear interpolation).
		Vec3<T> Interpolate(const Vec3<T> &ijk, int spatial_order) const;
    T InterpolateComponent(const Vec3<T> &ijk, int component,
                           int spatial_order = 1) const;

		// IO
		void Load(std::string file_name);
		void WriteVDB(std::string file_name) const;
		void Write(std::string file_name, bool text=false) const;
    void WriteCheckpoint() const {
      std::string checkpoint = "/tmp/" + name_ + "_" +
                               std::to_string(metagrid_.rank());
      Write(checkpoint, false);
    }

    // Send and receive data.
    void SendGhostData(
      canary::CanaryTaskContext *task_context, int ghost_width) const;
    int ReceiveGhostData(
      canary::CanaryTaskContext *task_context, int ghost_width);

    // Serialization and deserialization methods.
    // Serialize ghost regions and return rank for the neighboring partition.
    // If there is no neighbor in the given direction (boundary), return -1.
    int SerializeGhost(
      int di, int dj, int dk, int ghost_width,
      canary::CanaryOutputArchive &archive) const;
    int SerializeGhost(
      const Coord neighbor_dir, int ghost_width,
      canary::CanaryOutputArchive &archive) const;
    void DeserializeGhost(
      int ghost_width, canary::CanaryInputArchive &archive);

    // Serialization and deserialization methods using archive.
    template <class Archive>
    void save(Archive &ar) const {
      ar(metagrid_, staggered_, voxel_len_, name_);
      ar(background_[0], background_[1], background_[2]);
      CoordBBox bbox = data_->evalActiveVoxelBoundingBox();
      SerializeDense(bbox, ar);
    }
    template <class Archive>
    void load(Archive &ar) {
      ar(metagrid_, staggered_, voxel_len_, name_);
      ar(background_[0], background_[1], background_[2]);
      InitializePartition(metagrid_, staggered_, voxel_len_, background_, name_);
      DeserializeHelper(ar);
    }

    // Extrpolate non-normal components into boundary.
    void ExtrapolateToBoundary();

  private:
    // Helper serialization/deserialization methods.
    void SerializeDense(
      const CoordBBox box, canary::CanaryOutputArchive &archive) const;
    void DeserializeHelper(canary::CanaryInputArchive &archive);

};  // class VectorGrid

}  // namespace common

#endif  // COMMON_VECTOR_GRID_H
