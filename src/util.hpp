//
// Created by nickl on 12/27/18.
//

#ifndef MESH_SIMPL_UTIL_HPP
#define MESH_SIMPL_UTIL_HPP

#include <cmath>
#include "types.hpp"

namespace MeshSimpl {

vec3d &operator/=(vec3d &lfs, double rhs);

vec3d operator-(const vec3d &a, const vec3d &b);

namespace Internal {

// Dot product of two vectors
inline double dot(const vec3d &a, const vec3d &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Cross product of two vectors
inline vec3d cross(const vec3d &a, const vec3d &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

// Geometry midpoint of two positions
inline vec3d midpoint(const vec3d &a, const vec3d &b) {
  return {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
}

// |X|
inline double magnitude(const vec3d &x) {
  return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

inline order next(order ord) { return (ord + 1) % 3; }

inline order prev(order ord) { return (ord + 2) % 3; }

}  // namespace Internal

}  // namespace MeshSimpl

#endif  // MESH_SIMPL_UTIL_HPP
