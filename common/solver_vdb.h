#ifndef COMMON_SOLVER_VDB_H
#define COMMON_SOLVER_VDB_H

#include <openvdb/math/ConjGradient.h>
#include <openvdb/tools/Morphology.h>
#include <openvdb/tools/PoissonSolver.h>
#include <openvdb/tree/Tree.h>
#include <openvdb/Types.h>
#include <openvdb/util/NullInterrupter.h>

#include "common/definitions.h"

namespace common {

typedef int32_t VIndex;

namespace pcg = openvdb::math::pcg;

  typedef pcg::SparseStencilMatrix<double, 7> LaplacianMatrix;
typedef pcg::IncompleteCholeskyPreconditioner<openvdb::tools::poisson::LaplacianMatrix> IncompleteCholeskyPreconditioner;
typedef pcg::JacobiPreconditioner<openvdb::tools::poisson::LaplacianMatrix> JacobiPreconditioner;
typedef pcg::State State;

// Poisson solve, mostly copied as is, from openvdb/tools/PoissonSolver.h.
template<typename PreconditionerType, typename TreeType, typename BoundaryOp>
typename TreeType::Ptr
SolveWithBoundaryConditionsAndPreconditioner(
		const TreeType &input, const BoundaryOp &boundary_op,
		State &state, bool staggered) {

	typedef typename TreeType::ValueType TreeValueT;
	typedef typename LaplacianMatrix::ValueType VecValueT;
	typedef typename pcg::Vector<VecValueT> VectorT;
	typedef typename TreeType::template ValueConverter<VIndex>::Type VIdxTreeT;
	typedef typename TreeType::template ValueConverter<bool>::Type MaskTreeT;

	openvdb::util::NullInterrupter interrupter;

	// Create a mapping from active voxels of the input tree to elements of a
	// vector, and then populate a vector with input values.
	// When distributing, this mapping is for private + ghost regions, that is,
	// all data in the local tree, but the code will mostly stay the same,
	// except that we'll need to expand/reimplement below methods without
	// leafmanager's (tbb's) parallel implementation.
	typename VIdxTreeT::ConstPtr idx_tree = openvdb::tools::poisson::
		createIndexTree(input);
	typename VectorT::Ptr b = openvdb::tools::poisson::
		createVectorFromTree<VecValueT>(input, *idx_tree);

	// Assemble Laplacian matrix. When distributing, this the methods should
	// operate on private + ghost region, and need to be reimplemented without
	// tbb's foreach.
	typename MaskTreeT::Ptr interior_mask(
			new MaskTreeT(*idx_tree, false, openvdb::TopologyCopy()));
	openvdb::tools::erodeVoxels(*interior_mask, 1, openvdb::tools::NN_FACE);
	LaplacianMatrix::Ptr laplacian =
		openvdb::tools::poisson::createISLaplacianWithBoundaryConditions(
			*idx_tree, *interior_mask, boundary_op, *b, staggered);
	laplacian->scale(-1.0);
	b->scale(-1.0);
  // std::cout << "Laplacian matrix :" << std::endl;
  // std::cout << laplacian->str() << std::endl;

	// Solve the poisson equation.
	typename VectorT::Ptr x(
			new VectorT(b->size(), openvdb::zeroVal<VecValueT>()));
	typename pcg::Preconditioner<VecValueT>::Ptr precond(
			new PreconditionerType(*laplacian));
	if (!precond->isValid()) {
		precond.reset( new pcg::
				JacobiPreconditioner<LaplacianMatrix>(*laplacian));
	}
	// This needs to change when distributed, to a parallel implementation.
	pcg::State result = pcg::solve(
			*laplacian, *b, *x, *precond, interrupter, state);

  state.iterations = result.iterations;

	// Populate the output tree with values from the solution vector.
	return openvdb::tools::poisson::createTreeFromVector<TreeValueT>(
			*x, *idx_tree, openvdb::zeroVal<TreeValueT>());
}  // SolveWithBoundaryConditionsAndPreconditioner

}  // namespace common

#endif  // COMMON_SOLVER_VDB_H
