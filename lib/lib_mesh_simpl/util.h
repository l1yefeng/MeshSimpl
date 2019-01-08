//
// Created by nickl on 12/27/18.
//

#ifndef MESH_SIMPL_UTIL_H
#define MESH_SIMPL_UTIL_H

#include "types.h"
#include <array>
#include <vector>
#include <cassert>
#include <cmath>

namespace MeshSimpl
{

namespace Internal
{

template<typename Position>
inline double dot(const vec3d& a, const Position& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

inline vec3d cross(const vec3d& a, const vec3d& b)
{
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}

inline vec3d midpoint(const std::vector<double>& a, const std::vector<double>& b)
{
    return {(a[0]+b[0])/2, (a[1]+b[1])/2, (a[2]+b[2])/2};
}

inline double magnitude(const vec3d& x) { return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

template<typename PositionSrc, typename PositionDest>
inline void copy_vertex_position(const PositionSrc& src, PositionDest& dest)
{
    std::copy(src.begin(), src.begin()+3, dest.begin());
}

// Calculate quadric Q = (A, b, c) = (nn', dn, d*d)
inline Quadric make_quadric(const vec3d& normal, double d)
{
    return {
        /* A row 1 */ normal[0]*normal[0], normal[0]*normal[1], normal[0]*normal[2],
        /* A row 2 */                      normal[1]*normal[1], normal[1]*normal[2],
        /* A row 3 */                                           normal[2]*normal[2],
        /*    b    */ normal[0]*d, normal[1]*d, normal[2]*d,
        /*    c    */ d*d};
}

// Compute quadric error Q(v) = vAv + 2bv + c
template<typename Position>
inline double q_error(const Quadric& q, const Position& v)
{
    return dot({dot({q[0], q[1], q[2]}, v),
                dot({q[1], q[3], q[4]}, v),
                dot({q[2], q[4], q[5]}, v)}, v)+dot({q[6], q[7], q[8]}, v)*2+q[9];
}

} // namespace MeshSimpl::Internal

} // namespace MeshSimpl

MeshSimpl::Internal::Quadric&
operator+=(MeshSimpl::Internal::Quadric& lfs, const MeshSimpl::Internal::Quadric& rhs)
{
    for (unsigned int i = 0; i < 10; ++i)
        lfs[i] += rhs[i];
    return lfs;
}

MeshSimpl::Internal::Quadric&
operator*=(MeshSimpl::Internal::Quadric& lfs, double rhs)
{
    for (auto& x : lfs)
        x *= rhs;
    return lfs;
}

MeshSimpl::Internal::Quadric
operator+(const MeshSimpl::Internal::Quadric& a, const MeshSimpl::Internal::Quadric& b)
{
    MeshSimpl::Internal::Quadric q;
    for (unsigned int i = 0; i < 10; ++i)
        q[i] = a[i]+b[i];
    return q;
}

MeshSimpl::vec3d&
operator/=(MeshSimpl::vec3d& lfs, double rhs)
{
    for (auto& x : lfs)
        x /= rhs;
    return lfs;
}

MeshSimpl::vec3d
operator-(const std::vector<double>& a, const std::vector<double>& b)
{
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

MeshSimpl::vec3d
operator-(const MeshSimpl::vec3d& a, const std::vector<double>& b)
{
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

#endif // MESH_SIMPL_UTIL_H
