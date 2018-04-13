#include <functional>
#include <limits>

#include "common/parallel_utils.h"
#include "canary/canary.h"

namespace common {

template<typename T, typename Function>
std::enable_if<
  std::is_same<std::decay_t<std::result_of_t<Function(T, T)>>, T>::value, void>
Reduce(canary::CanaryApplication::VariableHandle<T> &local_handle,
       canary::CanaryApplication::VariableHandle<T> &global_handle,
       canary::CanaryApplication *app,
       Function combiner, T default_value) {
  // First scatter variables to rank 0.
  app->ReadAccess(local_handle);
  app->Scatter([=](canary::CanaryTaskContext *task_context) {
    task_context->Scatter(0, task_context->ReadVariable(local_handle));
  });

  // Then gather and reduce the value using the combiner.
  app->WriteAccess(global_handle);
  app->Gather([=](canary::CanaryTaskContext *task_context) -> int {
    int num_gathers = task_context->GetScatterParallelism();
    EXPECT_GATHER_SIZE(num_gathers);
    T *global = task_context->WriteVariable(global_handle);
    *global = default_value;
    *global = task_context->Reduce(*global, combiner);
    return 0;
  });

  // Finally broadcast reduced value.
  app->ReadAccess(global_handle);
  app->Scatter([=](canary::CanaryTaskContext *task_context) {
    task_context->Broadcast(task_context->ReadVariable(global_handle));
  });
  app->WriteAccess(local_handle);
  app->Gather([=](canary::CanaryTaskContext *task_context) -> int {
    EXPECT_GATHER_SIZE(1);
    task_context->GatherSingle(task_context->WriteVariable(local_handle));
    return 0;
  });
}  // Reduce

namespace {

template<typename T>
T sum(const T a, const T b) {
  return a+b;
}  // sum

template<typename T>
T max(const T a, const T b) {
  if (a > b) {
    return a;
  } else {
    return b;
  }
}  // max

template<typename T>
T min(const T a, const T b) {
  if (a <= b) {
    return a;
  } else {
    return b;
  }
}  // max

}  // anonymous

template<typename T>
void ReduceSum(canary::CanaryApplication::VariableHandle<T> &local_handle,
               canary::CanaryApplication::VariableHandle<T> &global_handle,
               canary::CanaryApplication *app) {
  Reduce(local_handle, global_handle, app, sum<T>, T(0));
}  // ReduceSum

template<typename T>
void ReduceMax(canary::CanaryApplication::VariableHandle<T> &local_handle,
               canary::CanaryApplication::VariableHandle<T> &global_handle,
               canary::CanaryApplication *app) {
  Reduce(local_handle, global_handle, app, max<T>,
         std::numeric_limits<T>::min());
}  // ReduceMax

template<typename T>
void ReduceMin(canary::CanaryApplication::VariableHandle<T> &local_handle,
               canary::CanaryApplication::VariableHandle<T> &global_handle,
               canary::CanaryApplication *app) {
  Reduce(local_handle, global_handle, app, min<T>,
         std::numeric_limits<T>::max());
}  // ReduceMax

template void ReduceSum(
  canary::CanaryApplication::VariableHandle<float> &local_handle,
  canary::CanaryApplication::VariableHandle<float> &global_handle,
  canary::CanaryApplication *app);
template void ReduceMax(
  canary::CanaryApplication::VariableHandle<float> &local_handle,
  canary::CanaryApplication::VariableHandle<float> &global_handle,
  canary::CanaryApplication *app);
template void ReduceMin(
  canary::CanaryApplication::VariableHandle<float> &local_handle,
  canary::CanaryApplication::VariableHandle<float> &global_handle,
  canary::CanaryApplication *app);

template void ReduceSum(
  canary::CanaryApplication::VariableHandle<double> &local_handle,
  canary::CanaryApplication::VariableHandle<double> &global_handle,
  canary::CanaryApplication *app);
template void ReduceMax(
  canary::CanaryApplication::VariableHandle<double> &local_handle,
  canary::CanaryApplication::VariableHandle<double> &global_handle,
  canary::CanaryApplication *app);
template void ReduceMin(
  canary::CanaryApplication::VariableHandle<double> &local_handle,
  canary::CanaryApplication::VariableHandle<double> &global_handle,
  canary::CanaryApplication *app);

}  // namespace common
