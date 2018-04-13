#include <string>
#include <sstream>

#include <openvdb/openvdb.h>

#include "common/primitives.h"

namespace common {

std::string ToString(Vec2i v) {
	std::stringstream s;
	s << v.x() << "," << v.y();
	return s.str();
}

std::string ToString(Vec2f v) {
	std::stringstream s;
	s << v.x() << "," << v.y();
	return s.str();
}

std::string ToString(Vec2d v) {
	std::stringstream s;
	s << v.x() << "," << v.y();
	return s.str();
}

std::string ToString(Vec3i v) {
	std::stringstream s;
	s << v.x() << "," << v.y() << "," << v.z();
	return s.str();
}

std::string ToString(Vec3f v) {
	std::stringstream s;
	s << v.x() << "," << v.y() << "," << v.z();
	return s.str();
}

std::string ToString(Vec3d v) {
	std::stringstream s;
	s << v.x() << "," << v.y() << "," << v.z();
	return s.str();
}

std::string ToString(Coord c) {
	std::stringstream s;
	s << c.x() << "," << c.y() << "," << c.z();
	return s.str();
}

std::string ToString(BBoxd b) {
	std::stringstream s;
	s << "{" << ToString(b.min()) << "},{" << ToString(b.max()) << "}";
	return s.str();
}

std::string ToString(CoordBBox b) {
	std::stringstream s;
	s << "{" << ToString(b.min()) << "},{" << ToString(b.max()) << "}";
	return s.str();
}

}  // namespace common
