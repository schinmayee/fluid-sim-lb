#include <cmath>
#include <iostream>

#include "common/solver_driver.h"

#include "canary/canary.h"
#include "common/parallel_utils.h"
#include "common/scalar_grid.h"
#include "common/solver.h"

namespace common {

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::Solve() {
  // First exchange ghost x_grid and ghost marker values.
  ExchangeGhostXAndMarker();

  // Assemble vectors and matrices, and build preconditioner.
  AssembleSolver();

  // Prepare pcg solver for the first iteration.
  PrepareSolver();

  // Global synchronization -- reduce rz and rr.
  Reduce_rr_rz();

  // Initialize current iteration number.
  InitializeIteration();

  // While loop.
  if (fixed_iterations_) {
    app_->Loop(max_iterations_);
  } else {
    WhileLoop();
  }

  //app_->TrackNeeded();

  // Exchange ghost p.
  ExchangeGhostP();

  // Local1.
  ComputeLocal1();

  // Global synchronization - reduce pAp.
  Reduce_pAp();

  // Local2.
  ComputeLocal2();

  // Global synchronization -- reduce rr and rz.
  Reduce_rr_rz();

  // Local3.
  ComputeLocal3();

  if (fixed_iterations_) {
    app_->EndLoop();
  } else {
    // Update number of iterations.
    NextIteration();
    // End while loop.
    app_->EndWhile();
  }
  PrintResidual();

  // Copy solution from vector to grid.
  CopySolution();
}

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::SetVariables(
  VariableHandle<ScalarGridT> *b_grid, VariableHandle<ScalarGridInt> *marker,
  VariableHandle<ScalarGridT> *x_grid) {
  b_grid_ = b_grid;
  marker_ = marker;
  x_grid_ = x_grid;
}  // SetVariables

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::DeclareVariables() {
  int num_partitions = partitions_[0]*partitions_[1]*partitions_[2];
  n_rows_ = new VariableHandle<int>(
    app_->DeclareVariable<int>(num_partitions));
  n_cols_ = new VariableHandle<int>(
    app_->DeclareVariable<int>(num_partitions));
  A_ = new VariableHandle<LaplacianMatrix>(
    app_->DeclareVariable<LaplacianMatrix>(num_partitions));
  idx_ = new VariableHandle<ScalarGridVIdx>(
    app_->DeclareVariable<ScalarGridVIdx>(num_partitions));
  b_ = new VariableHandle<VectorT>(
    app_->DeclareVariable<VectorT>(num_partitions));
  x_ = new VariableHandle<VectorT>(
    app_->DeclareVariable<VectorT>(num_partitions));
  r_ = new VariableHandle<VectorT>(
    app_->DeclareVariable<VectorT>(num_partitions));
  M_ = new VariableHandle<PreconditionerT>(
    app_->DeclareVariable<PreconditionerT>(num_partitions));
  p_ = new VariableHandle<VectorT>(
    app_->DeclareVariable<VectorT>(num_partitions));
  z_ = new VariableHandle<VectorT>(
    app_->DeclareVariable<VectorT>(num_partitions));
  Ap_ = new VariableHandle<VectorT>(
    app_->DeclareVariable<VectorT>(num_partitions));
  rz_local_ = new VariableHandle<T>(
    app_->DeclareVariable<T>(num_partitions));
  rz_old_local_ = new VariableHandle<T>(
    app_->DeclareVariable<T>(num_partitions));
  rr_local_ = new VariableHandle<T>(
    app_->DeclareVariable<T>(num_partitions));
  pAp_local_ = new VariableHandle<T>(
    app_->DeclareVariable<T>(num_partitions));
  rz_global_ = new VariableHandle<T>(
    app_->DeclareVariable<T>(1));
  rr_global_ = new VariableHandle<T>(
    app_->DeclareVariable<T>(1));
  pAp_global_ = new VariableHandle<T>(
    app_->DeclareVariable<T>(1));
  num_iterations_ = new VariableHandle<int>(
    app_->DeclareVariable<int>(1));
}  // DeclareVariables

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::AssignWorkerForGlobals(int idx) {
  std::vector<int> global = {idx};
  bool success = rz_global_->RecordInitialPartitionPlacement(app_, global);
  success &= rr_global_->RecordInitialPartitionPlacement(app_, global);
  assert(success);
  success &= pAp_global_->RecordInitialPartitionPlacement(app_, global);
  success &= num_iterations_->RecordInitialPartitionPlacement(app_, global);
  assert(success);
}  // AssignWorkerForGlobals

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::AssignWorkersForPartitions(
  const std::vector<int> &idxs) {
  bool success = n_rows_->RecordInitialPartitionPlacement(app_, idxs);
  success &= n_cols_->RecordInitialPartitionPlacement(app_, idxs);
  success &= b_grid_->RecordInitialPartitionPlacement(app_, idxs);
  success &= marker_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= x_grid_->RecordInitialPartitionPlacement(app_, idxs);
  success &= A_->RecordInitialPartitionPlacement(app_, idxs);
  success &= idx_->RecordInitialPartitionPlacement(app_, idxs);
  success &= b_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= x_->RecordInitialPartitionPlacement(app_, idxs);
  success &= r_->RecordInitialPartitionPlacement(app_, idxs);
  success &= M_->RecordInitialPartitionPlacement(app_, idxs);
  success &= p_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= z_->RecordInitialPartitionPlacement(app_, idxs);
  success &= Ap_->RecordInitialPartitionPlacement(app_, idxs);
  success &= rz_local_->RecordInitialPartitionPlacement(app_, idxs);
  success &= rz_old_local_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
  success &= rr_local_->RecordInitialPartitionPlacement(app_, idxs);
  success &= pAp_local_->RecordInitialPartitionPlacement(app_, idxs);
  assert(success);
}  // AssignWorkersForPartitions

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::SetMigratableVariables() {
  n_rows_->SetUpdatePlacement(app_, true);
  n_cols_->SetUpdatePlacement(app_, true);
  b_grid_->SetUpdatePlacement(app_, true);
  marker_->SetUpdatePlacement(app_, true);
  x_grid_->SetUpdatePlacement(app_, true);
  A_->SetUpdatePlacement(app_, true);
  idx_->SetUpdatePlacement(app_, true);
  b_->SetUpdatePlacement(app_, true);
  x_->SetUpdatePlacement(app_, true);
  r_->SetUpdatePlacement(app_, true);
  M_->SetUpdatePlacement(app_, true);
  p_->SetUpdatePlacement(app_, true);
  z_->SetUpdatePlacement(app_, true);
  Ap_->SetUpdatePlacement(app_, true);
  rz_local_->SetUpdatePlacement(app_, true);
  rz_old_local_->SetUpdatePlacement(app_, true);
  rr_local_->SetUpdatePlacement(app_, true);
  pAp_local_->SetUpdatePlacement(app_, true);
}  // SetMigratableVariables

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::ExchangeGhostXAndMarker() {
  // Send x_grid data.
  app_->ReadAccess(*x_grid_);
  app_->Scatter(
    [=](canary::CanaryTaskContext *task_context) {
      const ScalarGridT &x_grid = task_context->ReadVariable(*x_grid_);
      x_grid.SendGhostData(task_context, 1);
    }, name_+"_SendGhostX");

  // Receive x_grid data.
  app_->WriteAccess(*x_grid_);
  app_->Gather(
    [=](canary::CanaryTaskContext *task_context) -> int {
      ScalarGridT *x_grid = task_context->WriteVariable(*x_grid_);
      int ret = x_grid->ReceiveGhostData(task_context, 1);
      return ret;
    }, name_+"_ReceiveGhostX");
  // Send marker data.
  app_->ReadAccess(*marker_);
  app_->Scatter(
    [=](canary::CanaryTaskContext *task_context) {
      const ScalarGridInt &marker = task_context->ReadVariable(*marker_);
      marker.SendGhostData(task_context, 1);
    }, name_+"_SendGhostMarker");

  // Receive marker data.
  app_->WriteAccess(*marker_);
  app_->Gather(
    [=](canary::CanaryTaskContext *task_context) -> int {
      ScalarGridInt *marker = task_context->WriteVariable(*marker_);
      int ret = marker->ReceiveGhostData(task_context, 1);
      return ret;
    }, name_+"_ReceiveGhostMarker");
}  // ExchangeGhostXAndMarker

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::AssembleSolver() {
  app_->ReadAccess(*marker_);
  app_->ReadAccess(*b_grid_);
  app_->ReadAccess(*x_grid_);
  app_->WriteAccess(*idx_);
  app_->WriteAccess(*A_);
  app_->WriteAccess(*b_);
  app_->WriteAccess(*x_);
  app_->WriteAccess(*M_);
  app_->WriteAccess(*n_rows_);
  app_->WriteAccess(*n_cols_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << "AssembleSolver start ...";

      // Get variables.
      ScalarGridInt &marker = const_cast<ScalarGridInt&>(
        task_context->ReadVariable(*marker_));
      ScalarGridT &b_grid = const_cast<ScalarGridT&>(
        task_context->ReadVariable(*b_grid_));
      ScalarGridT &x_grid = const_cast<ScalarGridT&>(
        task_context->ReadVariable(*x_grid_));
      ScalarGridVIdx *idx = task_context->WriteVariable(*idx_);
      LaplacianMatrix *A = task_context->WriteVariable(*A_);
      VectorT *b = task_context->WriteVariable(*b_);
      VectorT *x = task_context->WriteVariable(*x_);
      PreconditionerT *M = task_context->WriteVariable(*M_);

      // Instantiate a solver.
      PCGSolver<ScalarGridT> solver;
      // Initialize solver metadata.
      solver.InitializePartition(marker.metagrid());
      // Set variables for assembly.
      solver.SetVariablesForAssembly(
        &marker, &b_grid, &x_grid, idx, A, b, x, M);
      // Assemble/initialize matrix, vectors, and preconditioner.
      solver.InitializeMatrixVectorsAndPreconditioner();
      // Get number of rows and columns.
      int *n_rows = task_context->WriteVariable(*n_rows_);
      int *n_cols = task_context->WriteVariable(*n_cols_);
      *n_rows = solver.n_rows();
      *n_cols = solver.n_cols();

      VLOG(1) << "AssembleSolver end ...";
    }, name_ + "_AssembleSolver");
}  // AssembleSolver

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::PrepareSolver() {
  app_->ReadAccess(*n_rows_);
  app_->ReadAccess(*n_cols_);
  app_->ReadAccess(*A_);
  app_->ReadAccess(*M_);
  app_->ReadAccess(*b_);
  app_->WriteAccess(*x_);
  app_->WriteAccess(*r_);
  app_->WriteAccess(*p_);
  app_->WriteAccess(*z_);
  app_->WriteAccess(*rr_local_);
  app_->WriteAccess(*rz_local_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << "PrepareSolver start ...";

      // Get variables.
      LaplacianMatrix &A = const_cast<LaplacianMatrix&>(
        task_context->ReadVariable(*A_));
      PreconditionerT &M = const_cast<PreconditionerT&>(
        task_context->ReadVariable(*M_));
      VectorT &b = const_cast<VectorT&>(
        task_context->ReadVariable(*b_));
      VectorT *x = task_context->WriteVariable(*x_);
      VectorT *r = task_context->WriteVariable(*r_);
      VectorT *p = task_context->WriteVariable(*p_);
      VectorT *z = task_context->WriteVariable(*z_);
      // Get scalars.
      T n_rows = task_context->ReadVariable(*n_rows_);
      T n_cols = task_context->ReadVariable(*n_cols_);
      T *rr = task_context->WriteVariable(*rr_local_);
      T *rz = task_context->WriteVariable(*rz_local_);
      // Instantiate a solver.
      PCGSolver<ScalarGridT> solver(n_rows, n_cols);
      // Set variables.
      solver.SetVariablesToPreparePCGSolver(&A, &M, &b, x, r, p, z);
      // Set scalars.
      solver.set_rr(*rr);
      solver.set_rz(*rz);
      // Invoke PreparePCGSolver.
      solver.PreparePCGSolver();
      // Save scalars.
      *rr = solver.rr();
      *rz = solver.rz();

      VLOG(1) << "PrepareSolver end ...";
    }, name_ + "_PrepareSolver");
}  // PrepareSolver

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::WhileLoop() {
  app_->ReadAccess(*num_iterations_);
  app_->ReadAccess(*rr_global_);
  app_->While(
    [=](canary::CanaryTaskContext *task_context) -> bool {
      //VLOG(1) << "While loop for PCG ...";
      int num_iterations = task_context->ReadVariable(*num_iterations_);
      T residual = task_context->ReadVariable(*rr_global_);
      assert(!std::isnan(residual));
      //VLOG(1) << "Solver iterations : " << num_iterations <<
      //             " and residual : " << residual;
      if (num_iterations >= max_iterations_ || residual < threshold_) {
        LOG(INFO) << "Solver done, after iterations = " << num_iterations <<
                     " residual = " << residual;
        return false;
      }
      return true;
    }, name_ + "_SolverWhileLoop");
}  // WhileLoop

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::Reduce_rr_rz() {
  ReduceSum<T>(*rr_local_, *rr_global_, app_);
  ReduceSum<T>(*rz_local_, *rz_global_, app_);
}  // Reduce_rr_zz

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::InitializeIteration() {
  app_->WriteAccess(*num_iterations_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << "InitializeIteration start ...";
      int *num_iterations = task_context->WriteVariable(*num_iterations_);
      *num_iterations = 0;
      VLOG(1) << "InitializeIteration end ...";
  }, name_ + "_SolverInitIterations");
}  // InitializeIteration

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::ExchangeGhostP() {
  // Send data.
  app_->ReadAccess(*p_);
  app_->ReadAccess(*idx_);
  app_->Scatter(
    [=](canary::CanaryTaskContext *task_context) {
      const VectorT &p = task_context->ReadVariable(*p_);
      const ScalarGridVIdx &idx = task_context->ReadVariable(*idx_);
      p.SendGhostData(task_context, 1, idx);
    });

  // Receive data.
  app_->WriteAccess(*p_);
  app_->ReadAccess(*idx_);
  app_->Gather(
    [=](canary::CanaryTaskContext *task_context) -> int {
      VectorT *p = task_context->WriteVariable(*p_);
      const ScalarGridVIdx &idx = task_context->ReadVariable(*idx_);
      int ret = p->ReceiveGhostData(task_context, 1, idx);
      return ret;
    });
}  // ExchangeGhostP

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::ComputeLocal1() {
  app_->ReadAccess(*n_rows_);
  app_->ReadAccess(*n_cols_);
  app_->ReadAccess(*A_);
  app_->ReadAccess(*p_);
  app_->WriteAccess(*Ap_);
  app_->WriteAccess(*pAp_local_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      //VLOG(1) << "ComputeLocal1 start ...";

      // Get variables.
      LaplacianMatrix &A = const_cast<LaplacianMatrix&>(
        task_context->ReadVariable(*A_));
      VectorT &p = const_cast<VectorT&>(task_context->ReadVariable(*p_));
      VectorT *Ap = task_context->WriteVariable(*Ap_);
      // Get scalars.
      T n_rows = task_context->ReadVariable(*n_rows_);
      T n_cols = task_context->ReadVariable(*n_cols_);
      T *pAp = task_context->WriteVariable(*pAp_local_);
      // Instantiate solver.
      PCGSolver<ScalarGridT> solver(n_rows, n_cols);
      // Set variables.
      solver.SetVariablesForComputeLocal1(&A, &p, Ap);
      // Set scalars.
      solver.set_pAp(*pAp);
      // Invoke ComputeLocal1.
      solver.ComputeLocal1();
      // Save scalars.
      *pAp = solver.pAp();

      //VLOG(1) << "ComputeLocal1 end ...";
    }, name_+"_ComputeLocal1");
}  // ComputeLocal1

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::Reduce_pAp() {
  ReduceSum<T>(*pAp_local_, *pAp_global_, app_);
}  // Reduce_pAp

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::ComputeLocal2() {
  app_->ReadAccess(*n_rows_);
  app_->ReadAccess(*n_cols_);
  app_->ReadAccess(*p_);
  app_->ReadAccess(*Ap_);
  app_->ReadAccess(*M_);
  app_->WriteAccess(*x_);
  app_->WriteAccess(*r_);
  app_->WriteAccess(*z_);
  app_->ReadAccess(*pAp_local_);
  app_->WriteAccess(*rz_local_);
  app_->WriteAccess(*rz_old_local_);
  app_->WriteAccess(*rr_local_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      //VLOG(1) << "ComputeLocal2 start ...";

      // Get variables.
      VectorT &p = const_cast<VectorT&>(task_context->ReadVariable(*p_));
      VectorT &Ap = const_cast<VectorT&>(task_context->ReadVariable(*Ap_));
      PreconditionerT &M =
        const_cast<PreconditionerT&>(task_context->ReadVariable(*M_));
      VectorT *x = task_context->WriteVariable(*x_);
      VectorT *r = task_context->WriteVariable(*r_);
      VectorT *z = task_context->WriteVariable(*z_);
      // Get scalars.
      T n_rows = task_context->ReadVariable(*n_rows_);
      T n_cols = task_context->ReadVariable(*n_cols_);
      T *rr = task_context->WriteVariable(*rr_local_);
      T *rz = task_context->WriteVariable(*rz_local_);
      T *rz_old = task_context->WriteVariable(*rz_old_local_);
      T pAp = task_context->ReadVariable(*pAp_local_);
      // Instantiate solver.
      PCGSolver<ScalarGridT> solver(n_rows, n_cols);
      // Set variables.
      solver.SetVariablesForComputeLocal2(x, &p, r, &Ap, &M, z);
      // Set scalars.
      solver.set_rr(*rr);
      solver.set_rz(*rz);
      solver.set_rz_old(*rz_old);
      solver.set_pAp(pAp);
      // Invoke ComputeLocal2.
      solver.ComputeLocal2();
      // Save scalars.
      *rr = solver.rr();
      *rz = solver.rz();
      *rz_old = solver.rz_old();

      //VLOG(1) << "ComputeLocal2 end ...";
    }, name_+"_ComputeLocal2");
}  // ComputeLocal2

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::ComputeLocal3() {
  app_->ReadAccess(*n_rows_);
  app_->ReadAccess(*n_cols_);
  app_->ReadAccess(*z_);
  app_->WriteAccess(*p_);
  app_->ReadAccess(*rz_local_);
  app_->ReadAccess(*rz_old_local_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      //VLOG(1) << "ComputeLocal3 start ...";

      // Get variables.
      VectorT &z = const_cast<VectorT&>(task_context->ReadVariable(*z_));
      VectorT *p = task_context->WriteVariable(*p_);
      // Get scalars.
      T n_rows = task_context->ReadVariable(*n_rows_);
      T n_cols = task_context->ReadVariable(*n_cols_);
      T rz = task_context->ReadVariable(*rz_local_);
      T rz_old = task_context->ReadVariable(*rz_old_local_);
      // Instantiate solver.
      PCGSolver<ScalarGridT> solver(n_rows, n_cols);
      // Set variables.
      solver.SetVariablesForComputeLocal3(p, &z);
      // Set scalars.
      solver.set_rz(rz);
      solver.set_rz_old(rz_old);
      // Invoke ComputeLocal3.
      solver.ComputeLocal3();
      // Save scalars -- none.
      
      //VLOG(1) << "ComputeLocal3 end ...";
    }, name_+"_ComputeLocal3");
}  // ComputeLocal3

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::PrintResidual() {
  app_->ReadAccess(*rr_global_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      LOG(INFO) << "Residual = " << task_context->ReadVariable(*rr_global_);
  }, "PrintResidual");  // Transform
}  // PrintResidual

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::NextIteration() {
  app_->WriteAccess(*num_iterations_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      int *num_iterations = task_context->WriteVariable(*num_iterations_);
      (*num_iterations)++;
    }, name_+"_NextIteration");
}  // NextIteration

template<class ScalarGridT>
void PCGDriver<ScalarGridT>::CopySolution() {
  app_->ReadAccess(*x_);
  app_->ReadAccess(*idx_);
  app_->ReadAccess(*n_rows_);
  app_->ReadAccess(*n_cols_);
  app_->WriteAccess(*x_grid_);
  app_->Transform(
    [=](canary::CanaryTaskContext *task_context) {
      VLOG(1) << "CopySolution start ...";
      // Get variables
      VectorT &x = const_cast<VectorT&>(task_context->ReadVariable(*x_));
      ScalarGridVIdx &idx = const_cast<ScalarGridVIdx&>(
        task_context->ReadVariable(*idx_));
      ScalarGridT *x_grid = task_context->WriteVariable(*x_grid_);
      // Get scalars.
      T n_rows = task_context->ReadVariable(*n_rows_);
      T n_cols = task_context->ReadVariable(*n_cols_);
      // Instantiate solver.
      PCGSolver<ScalarGridT> solver(n_rows, n_cols);
      // Initialize metadata.
      solver.InitializePartition(idx.metagrid());
      // Set variables.
      solver.SetVariablesToCopySolution(&x,  &idx, x_grid);
      // Invoke CopySolution.
      solver.CopySolution();
      VLOG(1) << "CopySolution end ...";
  });
}

template class PCGDriver< ScalarGrid<3, 3, float> >;
template class PCGDriver< ScalarGrid<3, 3, double> >;

}
