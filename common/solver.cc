#include <iostream>
#include <openvdb/math/ConjGradient.h>
#include <openvdb/tools/PoissonSolver.h>

#include "canary/canary.h"
#include "common/solver.h"
#include "common/scalar_grid.h"

#define LARGE_FLOAT 1E5

namespace common {

namespace pcg = openvdb::math::pcg;

template<class ScalarGridT>
void Solver<ScalarGridT>::InitializePartition(const MetaGrid &metagrid) {
  metagrid_valid_ = true;
  metagrid_ = metagrid;
}  // InitializePartition

template<class ScalarGridT>
void Solver<ScalarGridT>::ResetVariables() {
  if (b_grid_) {
    b_grid_ = nullptr;
  }
  if (marker_) {
    marker_ = nullptr;
  }
  if (x_grid_) {
    x_grid_ = nullptr;
  }
  if (A_) {
    A_ = nullptr;
  }
  if (idx_) {
    idx_ = nullptr;
  }
  if (x_) {
    x_ = nullptr;
  }
  if (r_) {
    r_ = nullptr;
  }
}  // ResetVariables

template<class ScalarGridT>
void Solver<ScalarGridT>::SetVariablesForAssembly(
  ScalarGridInt *marker, ScalarGridT *b_grid, ScalarGridT *x_grid,
  ScalarGridVIdx *idx, LaplacianMatrix *A, VectorT *b, VectorT *x,
  PreconditionerT *M) {
  marker_ = marker;
  b_grid_ = b_grid;
  x_grid_ = x_grid;
  idx_    = idx;
  A_      = A;
  b_      = b;
  x_      = x;
  M_      = M;
}  // SetVariablesForAssembly

template<class ScalarGridT>
void Solver<ScalarGridT>::SetVariablesToCopySolution(
  VectorT *x, ScalarGridVIdx *idx, ScalarGridT *x_grid) {
  x_ = x;
  idx_ = idx;
  x_grid_ = x_grid;
}  // SetVariablesToCopySolution

template<class ScalarGridT>
void Solver<ScalarGridT>::InitializeMatrixVectorsAndPreconditioner() {
  InitializeMatrixAndVectors();
  BuildPreconditioner();
}  // InitializeMatrixVectorsAndPreconditioner

template<class ScalarGridT>
void Solver<ScalarGridT>::InitializeMatrixAndVectors() {
  CHECK(metagrid_valid_);
  CoordBBox box_interior = metagrid_.bbox();
  CoordBBox box_total = box_interior.expandBy(1);
  n_rows_ = marker_->Count(box_interior, FLUID);
  n_cols_ = marker_->Count(box_total, FLUID);
  ScalarGridBool local_mask;
  local_mask.InitializePartition(metagrid_, 1.f, false);
  CreateIndexAndMaskTree(local_mask);
  b_->Resize(n_rows_);
  CreateVectorFromGrid(*b_grid_, box_interior, *b_);
  x_->Resize(n_cols_);
  CreateVectorFromGrid(*x_grid_, box_total, *x_);
  A_->Resize(n_rows_, n_cols_);
  CreateLaplacianMatrix(local_mask);
}  // InitializeMatrixAndVectors

template<class ScalarGridT>
void Solver<ScalarGridT>::CopySolution() {
  x_grid_->Clear();
  CHECK(metagrid_valid_);
  CoordBBox box_interior = metagrid_.bbox();
  CreateGridFromVector(*x_, box_interior, *x_grid_);
}  // CopySolution

template<class ScalarGridT>
void Solver<ScalarGridT>::BuildPreconditioner() {
  if (preconditioner_type_ == IC) {
    M_->Reset(typename PreconditionerT::PreconditionerT::Ptr(new pcg::
      IncompleteCholeskyPreconditioner<typename LaplacianMatrix::MatrixT>(
        *(A_->data()))));
  } else {
    M_->Reset(typename PreconditionerT::PreconditionerT::Ptr(new pcg::
      JacobiPreconditioner<typename LaplacianMatrix::MatrixT>(
        *(A_->data()))));
  }
}  // BuildPreconditioner

template<class ScalarGridT>
void Solver<ScalarGridT>::CreateIndexAndMaskTree(ScalarGridBool &local_mask) {
  CHECK(metagrid_valid_);
  // Create/initialize idx tree.
  idx_->InitializePartition(metagrid_);
  CoordBBox box_interior = metagrid_.bbox();
  CoordBBox box_total = box_interior.expandBy(1);
  const Coord start_interior = box_interior.min();
  const Coord end_interior   = box_interior.max();
  const Coord start_total = box_total.min();
  const Coord end_total   = box_total.max();
  int used_interior = 0;
  int used_exterior = n_rows_;
  for (int i = start_total[0]; i <= end_total[0]; ++i) {
    for (int j = start_total[1]; j <= end_total[1]; ++j) {
      for (int k = start_total[2]; k <= end_total[2]; ++k) {
        Coord ijk(i,j,k);
        if ( !((marker_->isOn(ijk)) && (marker_->get(ijk) == FLUID)) ) {
          continue;
        }
        if (i >= start_interior[0] && i <= end_interior[0] &&
            j >= start_interior[1] && j <= end_interior[1] &&
            k >= start_interior[2] && k <= end_interior[2]) {
          assert(used_interior < n_rows_);
          idx_->set(ijk, used_interior);
          local_mask.set(ijk, true); 
          used_interior++;
        } else {
          assert(used_exterior < n_cols_);
          idx_->set(ijk, used_exterior);
          used_exterior++;
        }
      }  // for k
    }  // for j
  }  // for i
  CHECK_EQ(used_interior, n_rows_);
  CHECK_EQ(used_exterior, n_cols_);
}  // CreateIndexAndMaskTree

template<class ScalarGridT>
void Solver<ScalarGridT>::CreateVectorFromGrid(
  const ScalarGridT &grid, const CoordBBox box, VectorT &vector) {
  auto data = vector.data();
  const Coord start = box.min();
  const Coord end = box.max();
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        Coord ijk(i,j,k);
        if (!idx_->isOn(ijk)) {
          continue;
        }
        VIndex index = idx_->get(ijk);
        CHECK(index < data->size());
        assert(index < data->size());
        data->at(index) = grid.get(ijk);
      }  // for k
    }  // for j
  }  // for i
}  // CreateVectorFromGrid

template<class ScalarGridT>
void Solver<ScalarGridT>::CreateGridFromVector(
  const VectorT &vector, const CoordBBox box, ScalarGridT &grid) {
  const typename VectorT::VectorT::Ptr data = vector.const_data();
  const Coord start = box.min();
  const Coord end = box.max();
  for (int i = start[0]; i <= end[0]; ++i) {
    for (int j = start[1]; j <= end[1]; ++j) {
      for (int k = start[2]; k <= end[2]; ++k) {
        Coord ijk(i,j,k);
        if (!idx_->isOn(ijk)) {
          continue;
        }
        VIndex index = idx_->get(ijk);
        CHECK(index < data->size());
        assert(index < data->size());
        grid.set(ijk, data->at(index));
      }  // for k
    }  // for j
  }  // for i
}  // CreateGridFromVector

template<class ScalarGridT>
void Solver<ScalarGridT>::CreateLaplacianMatrix(
  const ScalarGridBool &local_mask) {
  auto data = A_->data();
  typedef pcg::SizeType SizeType;
  const int num_offsets = 6;
  Coord offsets[num_offsets] = { Coord(-1,0,0), Coord(1,0,0),
                                 Coord(0,-1,0), Coord(0,1,0),
                                 Coord(0,0,-1), Coord(0,0,1) };
  for (typename ScalarGridBool::GridT::ValueOnCIter iter =
       local_mask.const_data()->cbeginValueOn(); iter.test(); ++iter) {
    const Coord ijk = iter.getCoord();
    const VIndex row_id = idx_->get(ijk);
    CHECK(row_id < n_rows_);
    CHECK_EQ(marker_->get(ijk), FLUID);
    assert(row_id < n_rows_);
    assert(marker_->get(ijk) == FLUID);
    T diagonal = 0;
    typename LaplacianMatrix::MatrixT::RowEditor row =
      data->getRowEditor(row_id);
    for (int n = 0; n < num_offsets; ++n) {
      const Coord nbr = ijk + offsets[n];
      const int nbr_marker = marker_->get(nbr);
      if (nbr_marker == FLUID || nbr_marker == AIR) {
        // Add 1 to diagonal.
        diagonal += 1;
      }
      if (nbr_marker == FLUID) {
        assert(idx_->isOn(nbr));
        // Set off-diagonal to -1.
        row.setValue(idx_->get(nbr), -1);
      }
    }  // for n < num_offsets
    // Set diagonal.
    row.setValue(row_id, diagonal);
  }  // for iter
}  // CreateLaplacianMatrix

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::ResetVariables() {
  Solver<ScalarGridT>::ResetVariables();
  if (p_) {
    p_ = nullptr;
  }
  if (z_) {
    z_ = nullptr;
  }
  if (Ap_) {
    Ap_ = nullptr;
  }
}  // ResetVariables

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::SetVariablesToPreparePCGSolver(
  LaplacianMatrix *A, PreconditionerT *M, VectorT *b, VectorT *x,
  VectorT *r, VectorT *p, VectorT *z) {
  Base::A_ = A;
  Base::b_ = b;
  Base::x_ = x;
  Base::r_ = r;
  Base::M_ = M;
  p_ = p;
  z_ = z;
}  // SetVariablesToPreparePCGSolver

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::PreparePCGSolver() {
  Base::r_->Resize(Base::n_rows_);
  z_->Resize(Base::n_rows_);
  p_->Resize(Base::n_cols_);
  auto A = Base::A_->data();
  auto x = Base::x_->data();
  auto b = Base::b_->data();
  auto M = Base::M_->data();
  auto r = Base::r_->data();
  auto z = z_->data();
  auto p = p_->data();
  // r = b - Ax
  pcg::internal::computeResidual(*A, *x, *b, *r);
  // z = M_inv r
  M->apply(*r, *z);
  // p = z (p and z have different sizes)
  for (int i = 0; i < Base::n_rows_; ++i) {
    p->at(i) = z->at(i);
  }
  // rz = r'z
  rz_ = r->dot(*z);
  // rr = r'r
  rr_ = r->dot(*r);
}  // PrepareCGSolver

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::SetVariablesForComputeLocal1(
  LaplacianMatrix *A, VectorT *p, VectorT *Ap) {
  Base::A_ = A;
  p_  = p;
  Ap_ = Ap;
}  // SetVariablesForComputeLocal1

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::ComputeLocal1() {
  Ap_->Resize(Base::n_rows_);
  auto A = Base::A_->data();
  auto p = p_->data();
  auto Ap = Ap_->data();
  // Ap = Ap
  A->vectorMultiply(*p, *Ap);
  pAp_ = 0;
  // pAp = p'Ap (p and AP have different sizes)
  for (int i = 0; i < Base::n_rows_; ++i) {
    pAp_ += p->at(i)*Ap->at(i);
  }
}  // ComputeLocal1

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::SetVariablesForComputeLocal2(
  VectorT *x, VectorT *p, VectorT *r, VectorT *Ap,
  PreconditionerT *M, VectorT *z) {
  Base::x_ = x;
  p_ = p;
  Base::r_ = r;
  Ap_ = Ap;
  Base::M_ = M;
  z_ = z;
}  // SetVariablesForComputeLocal2

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::ComputeLocal2() {
  auto &x  = *(Base::x_->data());
  auto &p  = *(p_->data());
  auto &r  = *(Base::r_->data());
  auto &Ap = *(Ap_->data());
  auto &M  = *(Base::M_->data());
  auto &z  = *(z_->data());
  // alpha = rz/pAp
  T alpha = rz_/pAp_;
  CHECK_LE(alpha, LARGE_FLOAT);
  // x = x + alpha * p
  pcg::internal::axpy(alpha, p, x, x);
  // r = r - alpha * Ap
  pcg::internal::axpy(-alpha, Ap, r, r);
  // z = M_inv * r
  M.apply(r, z);
  // rr = r r'r
  rr_ = r.dot(r);
  // rz_old = rz
  rz_old_ = rz_;
  // rznew = r'z
  rz_ = r.dot(z);
}  // ComputeLocal2

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::SetVariablesForComputeLocal3(
  VectorT *p, VectorT *z) {
  p_ = p;
  z_ = z;
}  // SetVariablesForComputeLocal3

template<class ScalarGridT>
void PCGSolver<ScalarGridT>::ComputeLocal3() {
  auto &p = *(p_->data());
  auto &z = *(z_->data());
  // beta = rz/rz_old
  T beta = rz_/rz_old_;
  CHECK_LE(beta, LARGE_FLOAT);
  // p = z + beta * p (p and z have different sizes)
  for (int i = 0; i < Base::n_rows_; ++i) {
    p[i] = z[i] + beta * p[i];
  }
}  // ComputeLocal3

template class Solver< ScalarGrid<3, 3, float> >;
template class Solver< ScalarGrid<3, 3, double> >;

template class PCGSolver< ScalarGrid<3, 3, float> >;
template class PCGSolver< ScalarGrid<3, 3, double> >;

}  // namespace common
