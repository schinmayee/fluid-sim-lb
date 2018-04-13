#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdint>
#include <numeric>
#include <vector>
#include "common/utils.h"

namespace common {

Timer::Timer() : elapsed_seconds_(0) {}

void Timer::Start() {
  start_ = std::chrono::system_clock::now();
}

void Timer::Pause() {
  end_ = std::chrono::system_clock::now();
  elapsed_seconds_ += (end_ - start_);
}

std::chrono::duration<double> Timer::Stop() {
	Pause();
  std::chrono::duration<double> result = elapsed_seconds_;
  Reset();
	return result;
}

void Timer::Reset() {
	elapsed_seconds_ = std::chrono::duration<double>(0.0);
}

std::chrono::duration<double> Timer::TimeElapsed() {
	return elapsed_seconds_;
}

// Implmentation taken from:
// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes.
template<typename T, bool ascending>
std::vector<size_t> SortIndexes(const std::vector<T> &v) {
  // Initialize original index locations.
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // Sort indexes based on comparing values in v, in descending order.
  if (ascending) {
    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) { return v[i1] < v[i2];});
  } else {
    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) { return v[i1] > v[i2];});
  }
  return idx;
}  // SortIndexes

template std::vector<size_t> SortIndexes<int, true>(const std::vector<int>&);
template std::vector<size_t> SortIndexes<float, true>(const std::vector<float>&);
template std::vector<size_t> SortIndexes<double, true>(const std::vector<double>&);
template std::vector<size_t> SortIndexes<int, false>(const std::vector<int>&);
template std::vector<size_t> SortIndexes<float, false>(const std::vector<float>&);
template std::vector<size_t> SortIndexes<double, false>(const std::vector<double>&);

}  // namespace common
