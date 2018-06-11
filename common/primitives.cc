#include <cmath>
#include <string>
#include <sstream>

#include <openvdb/openvdb.h>

#include "canary/canary.h"
#include "common/primitives.h"

namespace common {

template <typename T>
int Sphere<T>::ParseFrom(
    const std::vector<std::string> &tokens, const int start_idx) {
  int id = start_idx;
  id = ParseVec3<T>(tokens, id, &center_);
  radius_ = Parse<T>(tokens[id++]);
  return id;
}  // Sphere::ParseFrom

template <typename T>
void Sphere<T>::Scale(const common::Coord global_dims) {
  center_[0] *= T(global_dims[0]);
  center_[1] *= T(global_dims[1]);
  center_[2] *= T(global_dims[2]);
  radius_    *= T(global_dims[0]);
}  // Sphere::Scale

template <typename T>
std::string Sphere<T>::ToString() const {
  std::ostringstream ss;
  ss << "Center : " << common::ToString(center_);
  ss << ", Radius : " << radius_;
  return ss.str();
}  // Sphere::ToString

template <typename T>
CoordBBox Cube<T>::get_coord_bbox() const {
  const T tolerance = 1e-2;
  for (int d = 0; d < 3; ++d) {
    CHECK_LE(start_[d], tolerance + floor(start_[d]))
      << "Make sure that box for cube is scaled properly!";
    CHECK_LE(ceil(end_[d]), end_[d] + tolerance)
      << "Make sure that box for cube is scaled properly!";
  }
  Coord start(start_[0], start_[1], start_[2]);
  Coord end(end_[0]-1, end_[1]-1, end_[2]-1);
  for (int d = 0; d < 3; ++d) {
    if (start[d] > end[d]) {
      end[d] = start[d];
    }
  }
  return CoordBBox(start, end);
}  // get_coord_bbox

template <typename T>
bool Cube<T>::IsInside(const Vec3<T> &p) const {
  bool is_inside = true;
  for (int d = 0; d < 3; ++d) {
    is_inside &= (p[d] >= start_[d]);
    is_inside &= (p[d] <= end_[d]);
  }  // for d
  return is_inside;
}  // IsInside

template <typename T>
int Cube<T>::ParseFrom(
    const std::vector<std::string> &tokens, const int start_idx) {
  int id = start_idx;
  id = ParseVec3<T>(tokens, id, &start_);
  id = ParseVec3<T>(tokens, id, &end_);
  return id;
}  // Cube::ParseFrom

template <typename T>
void Cube<T>::Scale(const common::Coord global_dims) {
  start_[0] *= T(global_dims[0]);
  start_[1] *= T(global_dims[1]);
  start_[2] *= T(global_dims[2]);
  end_[0] *= T(global_dims[0]);
  end_[1] *= T(global_dims[1]);
  end_[2] *= T(global_dims[2]);
}  // Cube::Scale

template <typename T>
std::string Cube<T>::ToString() const {
  std::ostringstream ss;
  ss << "Start : " << common::ToString(start_);
  ss << ", End : " << common::ToString(end_);
  return ss.str();
}  // Sphere::ToString

template class Sphere<float>;
template class Sphere<double>;
template class Cube<float>;
template class Cube<double>;

}  // namespace common
