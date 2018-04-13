#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdint>
#include <vector>

namespace common {

#ifndef DEBUG

inline void Debug(FILE *fp, const char *format, ...) {}

#else  // NODEBUG

inline void Debug(FILE *fp, const char *format, ...) {
	va_list argptr;
	va_start(argptr, format);
	vfprintf(fp, format, argptr);
	va_end(argptr);
}

#endif  // DEBUG

class Timer {
	private:
    std::chrono::system_clock::time_point start_;
    std::chrono::system_clock::time_point end_;
    std::chrono::duration<double> elapsed_seconds_;
	public:
		Timer();
		void Start();
		void Pause();
		std::chrono::duration<double> Stop();
		void Reset();
		std::chrono::duration<double> TimeElapsed();
};

/*
 * Helper functions.
 */
template <typename T, bool ascending=false>
std::vector<size_t> SortIndexes(const std::vector<T> &v);

}  // namespace common

#endif  // COMMON_UTILS_H
