#include "common/matrix.h"

namespace common {

// Laplacian
template class Matrix<float, 7>;
template class Matrix<double, 7>;

// Triangular
template class Matrix<float, 4>;
template class Matrix<double, 4>;

// JacobiPreconditioner
template class JacobiPreconditioner<float>;
template class JacobiPreconditioner<double>;

// IncompleteCholeskyPreconditioner
template class IncompleteCholeskyPreconditioner<float>;
template class IncompleteCholeskyPreconditioner<double>;

}  // namespace common
