#ifndef COMMON_MATRIX_H
#define COMMON_MATRIX_H

#include <openvdb/math/ConjGradient.h>

#include "common/primitives.h"
#include "canary/canary.h"

namespace common {

template<typename T, int SIZE>
class Matrix {
  public:
    // Common types.
    typedef typename openvdb::math::pcg::SparseStencilMatrix<T, SIZE> MatrixT;
    typedef T ValueType;

  private:
    // Members.
    int num_rows_ = 0, num_cols_ = 0;
    typename MatrixT::Ptr data_;

  public:
    // Reset.
    inline void Reset() {
      num_rows_ = 0;
      num_cols_ = 0;
      data_.reset();
    }
    inline void Reset(typename MatrixT::Ptr data) {
      data_ = data;
    }
    inline void Resize(int num_rows, int num_cols) {
      Reset();
      num_rows_ = num_rows;
      num_cols_ = num_cols;
      data_ = typename MatrixT::Ptr(new MatrixT(num_rows, num_cols));
    }
    inline int num_rows() { return num_rows_; }

    // Accessors.
    inline typename MatrixT::Ptr data() { return data_; } 

    // Archive serialize.
    template<class Archive>
    void load(Archive &ar) {
      VLOG(1) << "Deserializing matrix ";
      int magic_number = 0;
      ar(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      ar(num_rows_, num_cols_);
      Resize(num_rows_, num_cols_);
      int num_elems;
      ar(num_elems);
      for (int i = 0; i < num_elems; ++i) {
        int r, c;
        T val;
        ar(r, c, val);
        data_->setValue(r, c, val);
      }
      magic_number = 0;
      ar(magic_number);
      CHECK_EQ(magic_number, kMagicNumber);
      VLOG(1) << "Deserialized matrix ";
    }
    template<class Archive>
    void save(Archive &ar) const {
      VLOG(1) << "Serializing matrix ";
      ar(kMagicNumber);
      ar(num_rows_, num_cols_);
      int num_elems = 0;
      for (int i = 0; i < num_rows_; ++i) {
        auto row = data_->getConstRow(i);
        for (auto iter = row.cbegin(); iter; ++iter) {
          num_elems++;
        }
      }
      ar(num_elems);
      int total = 0;
      for (int i = 0; i < num_rows_; ++i) {
        auto row = data_->getConstRow(i);
        for (auto iter = row.cbegin(); iter; ++iter) {
          int c = iter.column();
          T val = *iter;
          ar(i, c, val);
          total++;
        }
      }
      CHECK_EQ(total, num_elems);
      ar(kMagicNumber);
      VLOG(1) << "Serialized matrix ";
    }

};  // class Matrix

template<typename T>
class Preconditioner {
  public:
    // Common types.
    typedef typename
      openvdb::math::pcg::Preconditioner<T> PreconditionerT;
    typedef T ValueType;

  private:
    // Members.
    typename PreconditionerT::Ptr data_;

  public:
    // Reset.
    inline void Reset() { data_.reset(); }
    inline void Reset(typename PreconditionerT::Ptr data) { data_ = data; }

    // Accessors.
    inline typename PreconditionerT::Ptr data() { return data_; } 

    // Serialize to be implemented by child class.


    // Archive serialize.
    template<class Archive>
    void serialize(Archive &ar) {
      VLOG(1) << "*** Serialize method unimplemented for Preconditioner ***";
    }
};  // class Preconditioner

template<typename T>
class JacobiPreconditioner : Preconditioner<T> {
  public:
    // Common types.
    typedef typename
      openvdb::math::pcg::JacobiPreconditioner<T> JacobiT;
    typedef T ValueType;
};  // class JacobiPreconditioner

template<typename T>
class IncompleteCholeskyPreconditioner : Preconditioner<T> {
  public:
    // Common types.
    typedef typename openvdb::math::pcg::SparseStencilMatrix<T, 7> MatrixT;
    typedef typename
      openvdb::math::pcg::IncompleteCholeskyPreconditioner<MatrixT>
      IncompleteCholeskyT;
    typedef T ValueType;
};  // class IncompleteCholeskyPreconditioner

}  // namespace common

#endif // COMMON_MATRIX_H
