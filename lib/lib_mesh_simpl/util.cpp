//
// Created by nickl on 1/8/19.
//

#include "util.h"

namespace MeshSimpl {

Internal::Quadric& operator+=(Internal::Quadric& lfs, const Internal::Quadric& rhs) {
    for (unsigned int i = 0; i < 10; ++i)
        lfs[i] += rhs[i];
    return lfs;
}

Internal::Quadric& operator*=(Internal::Quadric& lfs, double rhs) {
    for (auto& x : lfs)
        x *= rhs;
    return lfs;
}

Internal::Quadric operator+(const Internal::Quadric& a, const Internal::Quadric& b) {
    Internal::Quadric q;
    for (unsigned int i = 0; i < 10; ++i)
        q[i] = a[i] + b[i];
    return q;
}

vec3d& operator/=(vec3d& lfs, double rhs) {
    for (auto& x : lfs)
        x /= rhs;
    return lfs;
}

vec3d operator-(const std::vector<double>& a, const std::vector<double>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

vec3d operator-(const vec3d& a, const std::vector<double>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

} // namespace MeshSimpl
