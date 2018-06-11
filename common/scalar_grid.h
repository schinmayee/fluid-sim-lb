#ifndef COMMON_SCALAR_GRID_H
#define COMMON_SCALAR_GRID_H

#include <cereal/archives/xml.hpp>
#include <iostream>
#include <string>

#include <openvdb/openvdb.h>

#include "canary/canary.h"
#include "common/metagrid.h"
#include "common/primitives.h"

namespace common {

template<Index N1, Index N2, typename T>
class ScalarGrid {
	public:
		typedef typename openvdb::tree::Tree3<T, N1, N2>::Type TreeT;
		typedef typename openvdb::Grid<TreeT> GridT;
    typedef T ValueType;
    template<typename ValueT>
    struct ValueConverter {
      using Type = ScalarGrid<N1, N2, ValueT>;
    };

		// constructor
		ScalarGrid();
		ScalarGrid(int x, int y, int z, float voxel_len = 1,
				       T background = T(0), std::string name = "data");
		ScalarGrid(Coord global_dims, float voxel_len = 1,
				       T background = T(0), std::string name = "data");

	private:
		// members
    MetaGrid metagrid_;
		T background_;  // background value
		float voxel_len_;  // length of each voxel
		std::string name_;
		typename GridT::Ptr data_;  // openvdb grid with data
		typename GridT::Accessor accessor_;  // grid accessor

	public:
		// methods

    void InitializePartition(
      const MetaGrid &metagrid,
      float voxel_len = 1, T background = T(0), std::string name = "data");

		inline int dims(size_t i) const { return metagrid_.dims(i); }
		inline float voxel_len() const { return voxel_len_; }
		inline T get(int i, int j, int k) const {
			return accessor_.getValue(Coord(i,j,k));
		}
		inline T get(Coord c) const {
			return accessor_.getValue(c);
		}
		inline void set(int i, int j, int k, T value) {
			accessor_.setValue(Coord(i,j,k), value);
		}
		inline void set(Coord c, T value) {
			accessor_.setValue(c, value);
		}
		inline bool isOn(int i, int j, int k) const {
			return accessor_.isValueOn(Coord(i,j,k));
		}
		inline bool isOn(Coord c) const {
			return accessor_.isValueOn(c);
		}
    inline void fill(const CoordBBox &box, T value, bool active=true) {
      accessor_.getTree()->denseFill(box, value, active);
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
		void Fill(T value);
		// Make the grid sparse by pruning.
		void Prune(T tolerance = T(0));
		// Copy to another grid.
		void CopyTo(ScalarGrid<N1, N2, T> &to);

    // Number of active cells within the given bounding box.
    int Count(CoordBBox box) const;
    // Number of cells equal to given value, in the provided bounding box.
    int Count(CoordBBox box, T value) const;

		// Interpolate to give posiition. Use spatial_order interpolation
		// (spatial_order = 1 means linear interpolation).
		T Interpolate(const Vec3<T> &ijk, int spatial_order) const;

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
      ar(metagrid_, background_, voxel_len_, name_);
      CoordBBox bbox = data_->evalActiveVoxelBoundingBox();
      Coord dims = bbox.dim();
      Coord local_dims = metagrid_.local_dims();
      SerializeDense(bbox, ar);
    }
    template <class Archive>
    void load(Archive &ar) {
      ar(metagrid_, background_, voxel_len_, name_);
      InitializePartition(metagrid_, voxel_len_, background_, name_);
      DeserializeHelper(ar);
    }

  private:
    // Helper serialization/deserialization methods.
    void SerializeDense(
      const CoordBBox box, canary::CanaryOutputArchive &archive) const;
    void DeserializeHelper(canary::CanaryInputArchive &archive);

};  // class ScalarGrid

}  // namespace common

#endif  // COMMON_SCALAR_GRID_H
