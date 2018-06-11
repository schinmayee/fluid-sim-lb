#ifndef PROJECTS_FLIP_COMMON_SOURCES_H
#define PROJECTS_FLIP_COMMON_SOURCES_H

#include <cereal/archives/xml.hpp>

#include "canary/canary.h"
#include "common/primitives.h"

namespace application {

template<typename T = float>
class Source {
  public:
    Source();
    explicit Source(
        const common::Cube<T> &box, const common::Vec3<T> &velocity);

    void set_box(const common::Cube<T> &box);
    void set_velocity(const common::Vec3<T> &velocity);
    void set_box(const common::Cube<T> &box1,
                 const common::Cube<T> &box2,
                 T period);
    void set_velocity(const common::Vec3<T> &velocity1,
                      const common::Vec3<T> &velocity2,
                      T period);
    void set_timed(T start, T end);

    common::Cube<T> get_box(T t = 0) const;
    common::Vec3<T> get_velocity(T t = 0) const;
    bool is_active(T t = 0) const;

    void ScaleBox(const common::Coord scale);
    void ScaleVelocity(T scale);

    std::string ToString() const;

  private:
    // Source bounding box --- the two regions between which the source is
    // oscillating.
    common::Cube<T> box1_, box2_;
    // Source velocity --- the two velocities between which the velocity is
    // oscillating.
    common::Vec3<T> velocity1_, velocity2_;
    // Oscillation periods for box and velocity.
    T box_period_ = 0;
    T velocity_period_ = 0;
    // Start and end time.
    T start_time_ = 0;
    T end_time_ = 0;
    // Is the source box oscillating?
    bool box_oscillating_ = false;
    // Is the source velocity oscillating?
    bool velocity_oscillating_ = false;
    // Is the source timed?
    bool timed_ = false;

  public:
    // Serialization/deserialization methods using archive.
    template<class Archive> void save(Archive &ar) const {
      const common::Vec3<T> &start1 = box1_.get_start();
      const common::Vec3<T> &end1   = box1_.get_end();
      const common::Vec3<T> &start2 = box2_.get_start();
      const common::Vec3<T> &end2   = box2_.get_end();
      ar(start1[0], start1[1], start1[2]);
      ar(end1[0], end1[1], end1[2]);
      ar(start2[0], start2[1], start2[2]);
      ar(end2[0], end2[1], end2[2]);
      ar(velocity1_[0], velocity1_[1], velocity1_[2]);
      ar(velocity2_[0], velocity2_[1], velocity2_[2]);
      ar(box_period_, velocity_period_, start_time_, end_time_);
      ar(box_oscillating_, velocity_oscillating_, timed_);
    }  // save
    template<class Archive> void load(Archive &ar) {
      common::Vec3<T> start1, end1, start2, end2;
      ar(start1[0], start1[1], start1[2]);
      ar(end1[0], end1[1], end1[2]);
      ar(start2[0], start2[1], start2[2]);
      ar(end2[0], end2[1], end2[2]);
      box1_.reset(start1, end1);
      box2_.reset(start2, end2);
      ar(velocity1_[0], velocity1_[1], velocity1_[2]);
      ar(velocity2_[0], velocity2_[1], velocity2_[2]);
      ar(box_period_, velocity_period_, start_time_, end_time_);
      ar(box_oscillating_, velocity_oscillating_, timed_);
    }  // load
};  // Source

}  // namespace application

#endif
