#ifndef COMMON_SOLVER_DRIVER_H
#define COMMON_SOLVER_DRIVER_H

#include <string>

#include "canary/canary.h"
#include "common/definitions.h"
#include "common/solver.h"

namespace common {

template<class ScalarGridT>
class PCGDriver {
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

    template<typename VariableType>
    using VariableHandle = canary::CanaryApplication::VariableHandle<VariableType>;

    PCGDriver() = delete;
    PCGDriver(canary::CanaryApplication *app, std::string name,
              int max_iterations, T threshold, bool fixed_iterations,
              PreconditionerType preconditioner_type,
              Coord partitions, Coord global_dims) :
      app_(app), name_(name),
      max_iterations_(max_iterations), threshold_(threshold),
      fixed_iterations_(fixed_iterations),
      preconditioner_type_(preconditioner_type),
      partitions_(partitions), global_dims_(global_dims), 
      n_rows_(nullptr), n_cols_(nullptr),
      b_grid_(nullptr), marker_(nullptr), x_grid_(nullptr), A_(nullptr),
      idx_(nullptr), b_(nullptr), x_(nullptr), r_(nullptr), M_(nullptr),
      p_(nullptr), z_(nullptr), Ap_(nullptr),
      rz_local_(nullptr), rz_old_local_(nullptr),
      rr_local_(nullptr), pAp_local_(nullptr),
      rz_global_(nullptr), rr_global_(nullptr), pAp_global_(nullptr),
      num_iterations_(nullptr) {}

    ~PCGDriver() {
      if (n_rows_)
        delete n_rows_;
      if (n_cols_)
        delete n_cols_;
      if (A_)
        delete A_;
      if (idx_)
        delete idx_;
      if (b_)
        delete b_;
      if (x_)
        delete x_;
      if (r_)
        delete r_;
      if (M_)
        delete M_;
      if (p_)
        delete p_;
      if (z_)
        delete z_;
      if (Ap_)
        delete Ap_;
      if (rz_local_)
        delete rz_local_;
      if (rz_old_local_)
        delete rz_old_local_;
      if (rr_local_)
        delete rr_local_;
      if (pAp_local_)
        delete pAp_local_;
      if (rz_global_)
        delete rz_global_;
      if (rr_global_)
        delete rr_global_;
      if (pAp_global_)
        delete pAp_global_;
      if (num_iterations_)
        delete num_iterations_;
    }  // ~PCGDriver

    // Set variables from application.
    void SetVariables(VariableHandle<ScalarGridT> *b_grid,
                      VariableHandle<ScalarGridInt> *marker,
                      VariableHandle<ScalarGridT> *x_grid);

    // Declare all canary variables here.
    void DeclareVariables();

    // Assign workers for partitions.
    void AssignWorkerForGlobals(int idx);
    void AssignWorkersForPartitions(const std::vector<int> &idxs);
    // Set which variables may be migrated.
    void SetMigratableVariables();

    // PCG computation -- invokes all transformations from building laplacian
    // matrix to solving using PCG.
    void Solve();
  
  private:
    canary::CanaryApplication *app_;
    std::string name_;
    int max_iterations_;
    T threshold_;
    bool fixed_iterations_;
    PreconditionerType preconditioner_type_;
    Coord partitions_;
    Coord global_dims_;

    // No need to serialize variable handles, they are for just defining the
    // program flow -- the DAG.
    VariableHandle<int> *n_rows_;
    VariableHandle<int> *n_cols_;

    VariableHandle<ScalarGridT> *b_grid_;
    VariableHandle<ScalarGridInt> *marker_;
    VariableHandle<ScalarGridT> *x_grid_;
    VariableHandle<LaplacianMatrix> *A_;
    VariableHandle<ScalarGridVIdx> *idx_;
    VariableHandle<VectorT> *b_;
    VariableHandle<VectorT> *x_;
    VariableHandle<VectorT> *r_;
    VariableHandle<PreconditionerT> *M_;
    VariableHandle<VectorT> *p_;
    VariableHandle<VectorT> *z_;
    VariableHandle<VectorT> *Ap_;

    // Scalars rz, rr, pAp are both local and global (since they are reduced).
    // Scalar rz_old is local only (saved from rz).
    VariableHandle<T> *rz_local_;
    VariableHandle<T> *rz_old_local_;
    VariableHandle<T> *rr_local_;
    VariableHandle<T> *pAp_local_;
    VariableHandle<T> *rz_global_;
    VariableHandle<T> *rr_global_;
    VariableHandle<T> *pAp_global_;
    VariableHandle<int> *num_iterations_;

    // Exchange ghost x and marker.
    void ExchangeGhostXAndMarker();

    // Assemble vectors and matrices and build block preconditioner.
    void AssembleSolver();

    // Prepare solver for the first iteration (compute all local quantities).
    void PrepareSolver();

    // While loop.
    void WhileLoop();

    // Initialize number of iterations.
    void InitializeIteration();

    // Reduce r'r and r'z.
    void Reduce_rr_rz();

    // Exchange ghost values for p.
    void ExchangeGhostP();

    // Compute local 1.
    void ComputeLocal1();

    // Reduce p'Ap.
    void Reduce_pAp();

    // Compute local 2.
    void ComputeLocal2();

    // Compute local 3.
    void ComputeLocal3();

    // Print residual.
    void PrintResidual();

    // Update iteration.
    void NextIteration();

    // Copy solution to grid.
    void CopySolution();
};

}  // namespace common

#endif  // COMMON_SOLVER_DRIVER_H
