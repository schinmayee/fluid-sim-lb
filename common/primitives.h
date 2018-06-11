/*
 * Defines primitives such as vectors, coordinates and coordinate boxes, and
 * some useful operations over those primitives.
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef COMMON_PRIMITIVES_H
#define COMMON_PRIMITIVES_H

#include <string>
#include <vector>

#include <openvdb/openvdb.h>

#include "canary/canary.h"

namespace common {

// OpenVDB already defines some useful arithmetic operations over these :)
// So just reuse.

// Values and co-ordinates.
typedef openvdb::Vec2i Vec2i;
typedef openvdb::Vec2f Vec2f;
typedef openvdb::Vec2d Vec2d;
typedef openvdb::Vec3i Vec3i;
typedef openvdb::Vec3f Vec3f;
typedef openvdb::Vec3d Vec3d;
typedef openvdb::Coord Coord;
typedef openvdb::BBoxd BBoxd;
typedef openvdb::CoordBBox CoordBBox;
typedef openvdb::Index Index;
typedef openvdb::Index64 Index64;
template <typename T>
using Vec2 = openvdb::math::Vec2<T>;
template <typename T>
using Vec3 = openvdb::math::Vec3<T>;

// Small matrices.
typedef openvdb::math::Mat3<float> Mat3f;
typedef openvdb::math::Mat3<double> Mat4f;
template <typename T>
using Mat3 = openvdb::math::Mat3<T>;

// Convert to string.

inline std::string ToString(const Vec2i &v) {
	std::stringstream s;
	s << v.x() << "," << v.y();
	return s.str();
}

inline std::string ToString(const Vec2f &v) {
	std::stringstream s;
	s << v.x() << "," << v.y();
	return s.str();
}

inline std::string ToString(const Vec2d &v) {
	std::stringstream s;
	s << v.x() << "," << v.y();
	return s.str();
}

inline std::string ToString(const Vec3i &v) {
	std::stringstream s;
	s << v.x() << "," << v.y() << "," << v.z();
	return s.str();
}

inline std::string ToString(const Vec3f &v) {
	std::stringstream s;
	s << v.x() << "," << v.y() << "," << v.z();
	return s.str();
}

inline std::string ToString(const Vec3d &v) {
	std::stringstream s;
	s << v.x() << "," << v.y() << "," << v.z();
	return s.str();
}

inline std::string ToString(const Coord &c) {
	std::stringstream s;
	s << c.x() << "," << c.y() << "," << c.z();
	return s.str();
}

inline std::string ToString(const BBoxd &b) {
	std::stringstream s;
	s << "{" << ToString(b.min()) << "},{" << ToString(b.max()) << "}";
	return s.str();
}

inline std::string ToString(const CoordBBox &b) {
	std::stringstream s;
	s << "{" << ToString(b.min()) << "},{" << ToString(b.max()) << "}";
	return s.str();
}

const int kMagicNumber = 47474747;

template <typename T>
inline T Parse(const std::string &token) {
  static_assert(sizeof(T) == -1, "Invalid call!");
  return T(0);
}

template<>
inline int Parse<int>(const std::string &token) {
  CHECK(token != "");
  return stoi(token);
}

template<>
inline float Parse<float>(const std::string &token) {
  CHECK(token != "");
  return stof(token);
}

template<>
inline double Parse<double>(const std::string &token) {
  CHECK(token != "");
  return stod(token);
}

// Parse Vec2, return next token id to parse.
template <typename T>
int ParseVec2(
    const std::vector<std::string> &tokens, int start_idx,
    Vec2<T> *result) {
  CHECK_LT(start_idx + 1, tokens.size());
  result->x() = Parse<T>(tokens[start_idx+0]);
  result->y() = Parse<T>(tokens[start_idx+1]);
  return start_idx+2;
}

// Parse Vec3, return next token id to parse.
template <typename T>
int ParseVec3(
    const std::vector<std::string> &tokens, int start_idx,
    Vec3<T> *result) {
  CHECK_LT(start_idx + 2, tokens.size());
  result->x() = Parse<T>(tokens[start_idx+0]);
  result->y() = Parse<T>(tokens[start_idx+1]);
  result->z() = Parse<T>(tokens[start_idx+2]);
  return start_idx+3;
}

template <typename T = float>
class Shape {
  public:
    // Parse from a set of tokens, return next token id to parse.
    virtual int ParseFrom(
        const std::vector<std::string> &tokens, const int start_idx) = 0;
    // Scale (multiply) values to change from normalized to unnormalized.
    virtual void Scale(const common::Coord global_dims) = 0;
    // Return a string descriptor.
    virtual std::string ToString() const = 0;

    // Inspect type.
    virtual bool IsSphere() const = 0;
    virtual bool IsCube() const = 0;
};  // class Shape

template <typename T = float>
class Sphere : public Shape<T> {
  private:
    Vec3<T> center_;
    T radius_;

  public:
    inline void set_center(const Vec3<T> &center) {
      center_ = center;
    }
    inline void set_radius(const T radius) {
      radius_ = radius;
    }
    inline Vec3<T> get_center() const {
      return center_;
    }
    inline T get_radius() const {
      return radius_;
    }

    virtual int ParseFrom(
        const std::vector<std::string> &tokens, const int start_idx);
    virtual void Scale(const common::Coord global_dims);
    virtual std::string ToString() const;

    virtual bool IsSphere() const {
      return true;
    }
    virtual bool IsCube() const {
      return false;
    }
};  // Sphere

template <typename T = float>
class Cube : public Shape<T> {
  private:
    common::Vec3<T> start_;
    common::Vec3<T> end_;

  public:
    Cube() {}
    Cube(const Vec3<T> &start, const Vec3<T> &end)
      : start_(start), end_(end) {}

    inline void reset(const Vec3<T> &start, const Vec3<T> &end) {
      start_ = start;
      end_ = end;
    }

    inline void set_start(const Vec3<T> &start) {
      start_ = start;
    }
    inline void set_end(const Vec3<T> &end) {
      end_ = end;
    }
    inline Vec3<T> get_start() const {
      return start_;
    }
    inline Vec3<T> get_end() const {
      return end_;
    }
    CoordBBox get_coord_bbox() const;
    bool IsInside(const Vec3<T> &p) const;

    virtual int ParseFrom(
        const std::vector<std::string> &tokens, const int start_idx);
    virtual void Scale(const common::Coord global_dims);
    virtual std::string ToString() const;

    virtual bool IsSphere() const {
      return false;
    }
    virtual bool IsCube() const {
      return true;
    }
};  // Cube

}  // namespace common

#endif  // COMMON_PRIMTIIVES_H
