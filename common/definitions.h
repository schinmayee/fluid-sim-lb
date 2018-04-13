#ifndef COMMON_DEFINITIONS_H
#define COMMON_DEFINITIONS_H

namespace common {

enum CellType { AIR, FLUID, SOLID, NUM_CELL_TYPES };
enum PreconditionerType { IC, IC0, JACOBI, NUM_PRECONDITIONER_TYPES };


}  // namespace common

#endif  // COMMON_DEFINITIONS_H
