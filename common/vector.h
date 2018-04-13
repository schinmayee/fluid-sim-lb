#ifndef COMMON_VECTOR_H
#define COMMON_VECTOR_H

#include <openvdb/math/ConjGradient.h>

#include "canary/canary.h"
#include "common/primitives.h"

namespace common {

template<typename T>
class Vector {
  public:
    // Common types.
    typedef typename openvdb::math::pcg::Vector<T> VectorT;
    typedef T ValueType;

  private:
    // Members.
    typename VectorT::Ptr data_;

  public:
    // Reset.
    inline void Reset() { data_.reset(); }
    inline void Reset(typename VectorT::Ptr data) { data_ = data; }

    // Intiialize to given size.
    inline void Resize(int n) {
      Reset();
      data_ = typename VectorT::Ptr(new VectorT(n));
    }
    inline void Resize(int n, T value) {
      Reset();
      data_ = typename VectorT::Ptr(new VectorT(n, value));
    }
    inline int size() const {
      if (data_)
        return data_->size();
      else
        return 0;
    }
    
    // Accessors.
    inline typename VectorT::Ptr data() { return data_; }
    inline typename VectorT::Ptr const_data() const { return data_; }

    // Send and receive data.
    template<class ScalarGridVIdx>
    void SendGhostData(
      canary::CanaryTaskContext *task_context, int ghost_width,
      const ScalarGridVIdx &idx) const;
    template<class ScalarGridVIdx>
    int ReceiveGhostData(
      canary::CanaryTaskContext *task_context, int ghost_width,
      const ScalarGridVIdx &idx);

    // Serialization and deserialization methods.
    // Serialize ghost regions and return rank for the neighboring partition.
    // If there is no neighbor in the given direction (boundary), return -1.
    template<class ScalarGridVIdx>
    int SerializeGhost(
      int di, int dj, int dk, int ghost_width, const ScalarGridVIdx &idx,
      canary::CanaryOutputArchive &archive) const;
    template<class ScalarGridVIdx>
    int SerializeGhost(
      const Coord neighbor_dir, int ghost_width, const ScalarGridVIdx &idx,
      canary::CanaryOutputArchive &archive) const;
    template<class ScalarGridVIdx>
    void DeserializeGhost(
      int ghost_width, const ScalarGridVIdx &idx,
      canary::CanaryInputArchive &archive);
    // Archive serialize.
    template<class Archive>
    void load(Archive &ar) {
      VLOG(1) << "Serializing vector ";
      int magic_number = 0;
      ar(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      int n;
      ar(n);
      Resize(n);
      for (int i = 0; i < n; ++i) {
        T val;
        ar(val);
        data_->at(i) = val;
      }
      magic_number = 0;
      ar(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      VLOG(1) << "Serialized vector ";
    }
    template<class Archive>
    void save(Archive &ar) const {
      VLOG(1) << "Deserializing vector ";
      ar(kMagicNumber);
      int n = size();
      ar(n);
      int total = 0;
      for (int i = 0; i < n; ++i) {
        ar(data_->at(i));
        total++;
      }
      CHECK_EQ(n, total);
      ar(kMagicNumber);
      VLOG(1) << "Deserialized vector ";
    }

  private:
    // Helper serialization/deserialization methods.
    template<class ScalarGridVIdx>
    void SerializeDense(
      const CoordBBox box, const ScalarGridVIdx &idx,
      canary::CanaryOutputArchive &archive) const;
    template<class ScalarGridVIdx>
    void DeserializeHelper(
      const ScalarGridVIdx &idx, canary::CanaryInputArchive &archive);

};  // class Vector

}  // namespace common

#endif  // COMMON_VECTOR_H
