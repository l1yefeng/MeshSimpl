//
// Created by nickl on 12/27/18.
//

#ifndef MESH_SIMPL_UTIL_H
#define MESH_SIMPL_UTIL_H

#include <array>
#include <vector>

namespace MeshSimpl
{

typedef unsigned int idx;
typedef std::array<double, 3> vec3d;
typedef std::array<unsigned int, 3> vec3u;
typedef std::vector<std::vector<double>> V;
typedef std::vector<std::vector<unsigned int>> F;

double dot(const vec3d& a, const vec3d& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
double dot(const vec3d& b, const std::vector<double>& a) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

namespace Internal
{

// A quadric Q consists of a symmetric 3x3 matrix A, a vec3 b, and a scalar c
typedef std::array<double, 10> Quadric;

// A convenient method converting std::vector to std::array (double * 3)
vec3d vec2arr(const std::vector<double>& v) { return {v[0], v[1], v[2]}; }

// Compute the optimal position: v = -inv(A)*b
vec3d optimal_v_pos(const Quadric& q)
{
    // computes the inverse of matrix A in quadric
    const double a_det = q[0]*(q[3]*q[5]-q[4]*q[4])
        -q[1]*(q[1]*q[5]-q[4]*q[2])
        +q[2]*(q[1]*q[4]-q[3]*q[2]);

    // FIXME: unsafe
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

MeshSimpl::Internal::Quadric
operator+(const MeshSimpl::Internal::Quadric& a, const MeshSimpl::Internal::Quadric& b)
{
    MeshSimpl::Internal::Quadric q(a);
    return q += b;
}

#endif // MESH_SIMPL_UTIL_H
