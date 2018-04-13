#ifndef COMMON_SOLVER_H
#define COMMON_SOLVER_H

#include <iostream>
#include <openvdb/math/ConjGradient.h>
#include <openvdb/tools/PoissonSolver.h>

#include "canary/canary.h"
#include "common/definitions.h"
#include "common/matrix.h"
#include "common/metagrid.h"
#include "common/primitives.h"
#include "common/vector.h"

namespace common {

template<class ScalarGridT>
class Solver {
  public:
    // Common types for use in solver.
    typedef int32_t VIndex;
    typedef typename ScalarGridT::ValueType T;
    typedef typename ScalarGridT::template ValueConverter<int32_t>::Type ScalarGridInt;
    typedef typename ScalarGridT::template ValueConverter<VIndex>::Type ScalarGridVIdx;
    typedef typename ScalarGridT::template ValueConverter<bool>::Type ScalarGridBool;
    typedef Vector<T> VectorT;
    typedef Matrix<T, 7> LaplacianMatrix;
    typedef Preconditioner<T> PreconditionerT;
    typedef IncompleteCholeskyPreconditioner<T> IncompleteCholeskyPreconditionerT;
    typedef JacobiPreconditioner<T> JacobiPreconditionerT;

    Solver() : preconditioner_type_(NUM_PRECONDITIONER_TYPES),
      b_grid_(nullptr), marker_(nullptr), x_grid_(nullptr), A_(nullptr),
      idx_(nullptr), x_(nullptr), r_(nullptr), M_(nullptr) {
      metagrid_valid_ = false;
      n_rows_ = 0;
      n_cols_ = 0;
    }
    Solver(int n_rows, int n_cols) : preconditioner_type_(NUM_PRECONDITIONER_TYPES),
      b_grid_(nullptr), marker_(nullptr), x_grid_(nullptr), A_(nullptr),
      idx_(nullptr), x_(nullptr), r_(nullptr), M_(nullptr) {
      metagrid_valid_ = false;
      n_rows_ = n_rows;
      n_cols_ = n_cols;
    }

    void InitializePartition(const MetaGrid &metagrid);
    
    template<class Archive> void serialize(Archive &ar) {
      std::cerr << "Solver base: serialize not implemented!" << std::endl;
      // TODO: Need to downcast and invoke methods on child class as
      // can't do virtual methods with template types.
    }  // serialize

    // Reset all variables in the solver.
    virtual void ResetVariables();

    // Set all variables needed for assembling matrix, vectors x and b,
    // and preconditioner.
    void SetVariablesForAssembly(
      ScalarGridInt *marker, ScalarGridT *b_grid, ScalarGridT *x_grid,
      ScalarGridVIdx *idx, LaplacianMatrix *A, VectorT *b, VectorT *x,
      PreconditionerT *M);

    // Set variables to copy solution.
    void SetVariablesToCopySolution(
      VectorT *x, ScalarGridVIdx *idx, ScalarGridT *x_grid);

    // Assemble vector, laplacian matrix, and build preconditioner.
    void InitializeMatrixVectorsAndPreconditioner();

    // Copy solution from vector to grid.
    void CopySolution();

  protected:
    // Stores region over which solver executes, and partition information.
    bool metagrid_valid_;
    MetaGrid metagrid_;
    // Preconditioner to use.
    PreconditionerType preconditioner_type_;
    // This is the rhs in Ax=b, in a mtrix form.
    ScalarGridT *b_grid_;
    // Marker used to assemble laplacian matrix.
    ScalarGridInt *marker_;
    // Value for x from previous step, in grid form.
    ScalarGridT *x_grid_;
    // A matrix.
    LaplacianMatrix *A_;
    // Idx grid storing mapping from grid to vector.
    ScalarGridVIdx *idx_;
    // Rhs in vector form.
    VectorT *b_;  // size N
    // Solution x in vector form.
    VectorT *x_;  // size N+G
    // Residual r in vector form.
    VectorT *r_;  // size N
    // Preconditioenr.
    PreconditionerT *M_;
    // Number of rows (internal cells) and columns (internal + ghost)
    int n_rows_;
    int n_cols_;

  private:

    // Initialize laplacian matrix, vectors, idx trees.
    void InitializeMatrixAndVectors();
    // Build preconditioner.
    void BuildPreconditioner();
    // Build index tree to map grid values to vectors, and a local mask which
    // is set wherever there is fluid in the interior.
    void CreateIndexAndMaskTree(ScalarGridBool &local_mask);
    // Linearize grid using index tree.
    void CreateVectorFromGrid(const ScalarGridT &grid, const CoordBBox box,
                              VectorT &vector);
    // Populate pressure grid with solution.
    void CreateGridFromVector(const VectorT &vector, const CoordBBox box,
                              ScalarGridT &grid);
    // Build Laplacian matrix.
    void CreateLaplacianMatrix(const ScalarGridBool &local_mask);

  public:
    int n_rows() { return n_rows_; }
    int n_cols() { return n_cols_; }
    void set_n_rows(int n_rows) { n_rows_ = n_rows; }
    void set_n_cols(int n_cols) { n_cols_ = n_cols; }
};  // class Solver

template<class ScalarGridT>
class PCGSolver : public Solver<ScalarGridT> {
  public:
    using Base = Solver<ScalarGridT>;
    typedef typename Base::VIndex VIndex;
    typedef typename Base::T T;
    typedef typename Base::ScalarGridInt ScalarGridInt;
    typedef typename Base::ScalarGridVIdx ScalarGridVIdx;
    typedef typename Base::ScalarGridBool ScalarGridBool;
    typedef typename Base::VectorT VectorT;
    typedef typename Base::LaplacianMatrix LaplacianMatrix;
    typedef typename Base::PreconditionerT PreconditionerT;
    typedef typename Base::IncompleteCholeskyPreconditionerT IncompleteCholeskyPreconditioenrT;
    typedef typename Base::JacobiPreconditionerT JacobiPreconditionerT;

    PCGSolver() : Solver<ScalarGridT>(), p_(nullptr), z_(nullptr) {}
    PCGSolver(int n_rows, int n_cols) : Solver<ScalarGridT>(n_rows, n_cols),
      p_(nullptr), z_(nullptr) {}

    // Reset all variables.
    virtual void ResetVariables();

    // Set all solve variables.
    void SetVariablesToPreparePCGSolver(
      LaplacianMatrix *A, PreconditionerT *M, VectorT *b, VectorT *x,
      VectorT *r, VectorT *p, VectorT *z);

    // Set variables for 1st iterative local stage.
    void SetVariablesForComputeLocal1(
      LaplacianMatrix *A, VectorT *p, VectorT *Ap);

    // Set variables for 2nd iterative local stage.
    void SetVariablesForComputeLocal2(
      VectorT *x, VectorT *p, VectorT *r, VectorT *Ap,
      PreconditionerT *M, VectorT *z);

    // Set variables for 3rd iterative local stage.
    void SetVariablesForComputeLocal3(VectorT *p, VectorT *z);

    // Prepare PCG solver -- run/intiailzie PCG variables before the iterations
    // start.
    void PreparePCGSolver();
    // First stage of iterative local computation.
    void ComputeLocal1();
    // Second stage of iterative local computation.
    void ComputeLocal2();
    // Third stage of iterative local computation.
    void ComputeLocal3();

    // Run PCG solve locally -- this does not work multinode. Use for testing
    // implementation on single node.
    void Solve();

  protected:
    // Direction p in vector form.
    VectorT *p_;  // size N+G
    // Conditioned residual z in vector form.
    VectorT *z_;  // size N
    // Product of A and p.
    VectorT *Ap_;  // size N
    // Dot product of r and z.
    T rz_;
    T rz_old_;
    // Dot product of r and r -- L2 norm squared.
    T rr_;
    // Product p'Ap.
    T pAp_;

  public:
    // Accessors
    T rz() { return rz_; }
    T rz_old() { return rz_old_; }
    T rr() { return rr_; }
    T pAp() { return pAp_; }
    void set_rz(T rz) { rz_ = rz; }
    void set_rz_old(T rz_old) { rz_old_ = rz_old; }
    void set_rr(T rr) { rr_ = rr; }
    void set_pAp(T pAp) { pAp_ = pAp; }
};  // class PCGSolver

}  // namespace common

#endif  // COMMON_SOLVER_H
