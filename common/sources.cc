#include <sstream>

#include "common/sources.h"

#include "canary/canary.h"
#include "common/primitives.h"

namespace application {

template<typename T>
Source<T>::Source() {
  box_oscillating_ = false;
  velocity_oscillating_ = false;
}  // Source

template<typename T>
Source<T>::Source(
    const common::Cube<T> &box, const common::Vec3<T> &velocity) {
  box1_ = box;
  velocity1_ = velocity;
  box_oscillating_ = false;
  velocity_oscillating_ = false;
}  // Source

template<typename T>
void Source<T>::set_box(const common::Cube<T> &box) {
  box1_ = box;
  box_oscillating_ = false;
}  // set_box

template<typename T>
void Source<T>::set_velocity(const common::Vec3<T> &velocity) {
  velocity1_ = velocity;
  velocity_oscillating_ = false;
}  // set_velocity

template<typename T>
void Source<T>::set_box(const common::Cube<T> &box1,
                        const common::Cube<T> &box2,
                        T period) {
  box1_ = box1;
  box2_ = box2;
  box_period_ = period;
  box_oscillating_ = true;
}  // set_box

template<typename T>
void Source<T>::set_velocity(const common::Vec3<T> &velocity1,
                             const common::Vec3<T> &velocity2,
                             T period) {
  velocity1_ = velocity1;
  velocity2_ = velocity2;
  velocity_period_ = period;
  velocity_oscillating_ = true;
}  // set_velocity

template<typename T>
void Source<T>::set_timed(T start, T end) {
  start_time_ = start;
  end_time_ = end;
  timed_ = true;
}  // set_timed

template<typename T>
common::Cube<T> Source<T>::get_box(T t) const {
  common::Cube<T> box = box1_;
  if (box_oscillating_) {
    common::Cube<T> box;
    common::Vec3<T> start, end;
    T w = M_PI * t / box_period_;
    for (int d = 0; d < 3; ++d) {
      T s1 = box1_.get_start()[d];
      T s2 = box2_.get_start()[d];
      T e1 = box1_.get_end()[d];
      T e2 = box2_.get_end()[d];
      start[d] = s1 + 0.5 * (s2-s1) * (cos(w)+1);
      end[d]   = e1 + 0.5 * (e2-e1) * (cos(w)+1);
    }  // for d
    box = common::Cube<T>(start, end);
  }
  box.set_end(box.get_end()-common::Vec3<T>(1));
  return box;
}  // get_box

template<typename T>
common::Vec3<T> Source<T>::get_velocity(T t) const {
  if (!velocity_oscillating_) {
    return velocity1_;
  }
  common::Vec3<T> velocity;
  T w = M_PI * t / box_period_;
  for (int d = 0; d < 3; ++d) {
    velocity[d] = velocity1_[d] +
                  0.5 * (velocity2_[d] - velocity1_[d]) * (cos(w) + 1);
  }  // for d
  return velocity;
}  // get_velocity

template<typename T>
bool Source<T>::is_active(T t) const {
  if (!timed_) {
    return true;
  }
  return (t >= start_time_ && t <= end_time_);
}  // is_active

template<typename T>
void Source<T>::ScaleBox(const common::Coord scale) {
  box1_.Scale(scale);
  box2_.Scale(scale);
}  // ScaleBox

template<typename T>
void Source<T>::ScaleVelocity(T scale) {
  velocity1_ *= scale;
  velocity2_ *= scale;
}  // ScaleVelocity

template<typename T>
std::string Source<T>::ToString() const {
  std::ostringstream ss;
  ss << "Dimensions: ";
  ss << box1_.ToString();
  if (box_oscillating_) {
    ss << " to ";
    ss << box2_.ToString();
  }
  ss << ", Velocity: ";
  ss << common::ToString(velocity1_);
  if (velocity_oscillating_) {
    ss << " to ";
    ss << common::ToString(velocity2_);
  }
  if (timed_) {
    ss << ", Active time: " << start_time_ << " to " << end_time_;
  }
  return ss.str();
}  // ToString

template class Source<float>;

}  // namespace application
