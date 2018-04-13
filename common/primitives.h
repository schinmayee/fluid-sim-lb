/*
 * Defines primitives such as vectors, coordinates and coordinate boxes, and
 * some useful operations over those primitives.
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef COMMON_PRIMITIVES_H
#define COMMON_PRIMITIVES_H

#include <string>

#include <openvdb/openvdb.h>

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
std::string ToString(Vec2i v);
std::string ToString(Vec2f v);
std::string ToString(Vec2d v);
std::string ToString(Vec3i v);
std::string ToString(Vec3f v);
std::string ToString(Vec3d v);
std::string ToString(Coord c);
std::string ToString(BBoxd b);
std::string ToString(CoordBBox b);

const int kMagicNumber = 47474747;

}  // namespace common

#endif  // COMMON_PRIMTIIVES_H
