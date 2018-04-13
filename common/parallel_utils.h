#ifndef PARALLEL_UTILS_H
#define PARALLEL_UTILS_H

#include "canary/canary.h"

namespace common {

template<typename T, typename Function>
std::enable_if<
  std::is_same<std::decay_t<std::result_of_t<Function(T, T)>>, T>::value, void>
Reduce(canary::CanaryApplication::VariableHandle<T> &local,
       canary::CanaryApplication::VariableHandle<T> &global,
       canary::CanaryApplication *app,
       Function combiner, T default_value);

template<typename T>
void ReduceSum(canary::CanaryApplication::VariableHandle<T> &local,
               canary::CanaryApplication::VariableHandle<T> &global,
               canary::CanaryApplication *app);

template<typename T>
void ReduceMax(canary::CanaryApplication::VariableHandle<T> &local,
               canary::CanaryApplication::VariableHandle<T> &global,
               canary::CanaryApplication *app);

template<typename T>
void ReduceMin(canary::CanaryApplication::VariableHandle<T> &local,
               canary::CanaryApplication::VariableHandle<T> &global,
               canary::CanaryApplication *app);

}  // namespace common

#endif  // PARALLEL_UTILS_H
