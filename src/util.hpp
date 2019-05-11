//
// Created by nickl on 12/27/18.
//

#ifndef LIB_MESH_SIMPL_UTIL_HPP
#define LIB_MESH_SIMPL_UTIL_HPP

#include <algorithm>
#include <cmath>
#include "neighbor.hpp"
#include "types.hpp"

namespace MeshSimpl {

Internal::Quadric& operator+=(Internal::Quadric& lfs,
                              const Internal::Quadric& rhs);

Internal::Quadric& operator*=(Internal::Quadric& lfs, double rhs);

Internal::Quadric operator+(const Internal::Quadric& a,
                            const Internal::Quadric& b);

vec3d& operator/=(vec3d& lfs, double rhs);

vec3d operator-(const vec3d& a, const vec3d& b);

namespace Internal {

// Dot product of two vectors
inline double dot(const vec3d& a, const vec3d& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Cross product of two vectors
inline vec3d cross(const vec3d& a, const vec3d& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

// Geometry midpoint of two positions
inline vec3d midpoint(const vec3d& a, const vec3d& b) {
  return {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
}

// |X|
inline double magnitude(const vec3d& x) {
  return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

// Calculate quadric Q = (A, b, c) = (nn', dn, d*d)
inline Quadric make_quadric(const vec3d& normal, double d) {
  return {normal[0] * normal[0], normal[0] * normal[1],
          normal[0] * normal[2], normal[1] * normal[1],
          normal[1] * normal[2], normal[2] * normal[2],
          normal[0] * d,         normal[1] * d,
          normal[2] * d,         d * d};
}

// Compute quadric error Q(v) = vAv + 2bv + c
inline double q_error(const Quadric& q, const vec3d& v) {
  return dot({dot({q[0], q[1], q[2]}, v), dot({q[1], q[3], q[4]}, v),
              dot({q[2], q[4], q[5]}, v)},
             v) +
         dot({q[6], q[7], q[8]}, v) * 2 + q[9];
}

bool sort_and_find_intersection(std::vector<idx>& values0,
                                std::vector<idx>& values1);

}  // namespace Internal

}  // namespace MeshSimpl

#endif  // LIB_MESH_SIMPL_UTIL_HPP
