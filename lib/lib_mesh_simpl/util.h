//
// Created by nickl on 12/27/18.
//

#ifndef MESH_SIMPL_UTIL_H
#define MESH_SIMPL_UTIL_H

#include <array>
#include <vector>
#include <cassert>
#include <cmath>

namespace MeshSimpl
{

typedef unsigned int idx;
typedef std::array<double, 3> vec3d;               // double
typedef std::array<idx, 3> vec3i;                  // idx
typedef std::array<idx, 2> vec2i;                  // idx
typedef std::vector<std::vector<double>> V;        // input/output vertex positions
typedef std::vector<std::vector<unsigned int>> F;  // input/output face indices

namespace Internal
{

// A quadric Q consists of a symmetric 3x3 matrix A, a vec3 b, and a scalar c
typedef std::array<double, 10> Quadric;

double dot(const vec3d& a, const vec3d& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

double dot(const vec3d& b, const std::vector<double>& a) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

vec3d cross(const vec3d& a, const vec3d& b)
{
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}

double magnitude(const vec3d& x) { return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

// Calculate quadric Q = (A, b, c) = (nn', dn, d*d)
Quadric make_quadric(const vec3d& normal, double d)
{
    return {
        /* A row 1 */ normal[0]*normal[0], normal[0]*normal[1], normal[0]*normal[2],
        /* A row 2 */                      normal[1]*normal[1], normal[1]*normal[2],
        /* A row 3 */                                           normal[2]*normal[2],
        /*    b    */ normal[0]*d, normal[1]*d, normal[2]*d,
        /*    c    */ d*d};
}

// Compute the optimal position: v = -inv(A)*b
vec3d optimal_v_pos(const Quadric& q)
{
    // computes the inverse of matrix A in quadric
    const double a_det = q[0]*(q[3]*q[5]-q[4]*q[4])
        -q[1]*(q[1]*q[5]-q[4]*q[2])
        +q[2]*(q[1]*q[4]-q[3]*q[2]);

    // FIXME: unsafe
    assert(a_det != 0);
    const double a_det_inv = 1.0/a_det;
    const std::array<double, 6> a_inv{(q[3]*q[5]-q[4]*q[4])*a_det_inv,
                                      (q[2]*q[4]-q[1]*q[5])*a_det_inv,
                                      (q[1]*q[4]-q[2]*q[3])*a_det_inv,
                                      (q[0]*q[5]-q[2]*q[2])*a_det_inv,
                                      (q[1]*q[2]-q[0]*q[4])*a_det_inv,
                                      (q[0]*q[3]-q[1]*q[1])*a_det_inv};

    const vec3d b{q[6], q[7], q[8]};
    return {-dot({a_inv[0], a_inv[1], a_inv[2]}, b),
            -dot({a_inv[1], a_inv[3], a_inv[4]}, b),
            -dot({a_inv[2], a_inv[4], a_inv[5]}, b)};
}

// Compute quadric error Q(v) = vAv + 2bv + c
double q_error(const Quadric& q, const vec3d& v, bool optimal)
{
    const vec3d b{q[6], q[7], q[8]};
    const double c = q[9];
    if (optimal)
        return dot(b, v)+c;

    return dot({dot({q[0], q[1], q[2]}, v),
                dot({q[1], q[3], q[4]}, v),
                dot({q[2], q[4], q[5]}, v)}, v)+dot(b, v)*2+c;
}

} // namespace MeshSimpl::Internal

} // namespace MeshSimpl

MeshSimpl::Internal::Quadric&
operator+=(MeshSimpl::Internal::Quadric& lfs, const MeshSimpl::Internal::Quadric& rhs)
{
    for (unsigned int i = 0; i < rhs.size(); ++i)
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
    MeshSimpl::Internal::Quadric q(a);
    return q += b;
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
