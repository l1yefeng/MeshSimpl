//
// Created by nickl on 12/27/18.
//

#ifndef LIB_MESH_SIMPL_UTIL_H
#define LIB_MESH_SIMPL_UTIL_H

#include "neighbor.hpp"
#include "types.hpp"
#include <cassert>
#include <cmath>

namespace MeshSimpl {

namespace Internal {

// Dot product of two vectors
template <typename Position> inline double dot(const vec3d& a, const Position& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Cross product of two vectors
inline vec3d cross(const vec3d& a, const vec3d& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

// Geometry midpoint of two positions
inline vec3d midpoint(const std::vector<double>& a, const std::vector<double>& b) {
    return {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
}

// |X|
inline double magnitude(const vec3d& x) {
    return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

template <typename PositionSrc, typename PositionDest>
inline void copy_vec3(const PositionSrc& src, PositionDest& dest) {
    std::copy(src.begin(), src.begin() + 3, dest.begin());
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
template <typename Position> inline double q_error(const Quadric& q, const Position& v) {
    return dot({dot({q[0], q[1], q[2]}, v), dot({q[1], q[3], q[4]}, v),
                dot({q[2], q[4], q[5]}, v)},
               v) +
           dot({q[6], q[7], q[8]}, v) * 2 + q[9];
}

inline order v_in_face(const F& indices, idx f, idx v) {
    assert(indices[f][0] == v || indices[f][1] == v || indices[f][2] == v);
    return indices[f][0] == v ? 0 : (indices[f][1] == v ? 1 : 2);
}

} // namespace Internal

Internal::Quadric& operator+=(Internal::Quadric& lfs, const Internal::Quadric& rhs);

Internal::Quadric& operator*=(Internal::Quadric& lfs, double rhs);

Internal::Quadric operator+(const Internal::Quadric& a, const Internal::Quadric& b);

vec3d& operator/=(vec3d& lfs, double rhs);

vec3d operator-(const std::vector<double>& a, const std::vector<double>& b);

vec3d operator-(const vec3d& a, const std::vector<double>& b);

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_UTIL_H
