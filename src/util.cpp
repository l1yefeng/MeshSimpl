//
// Created by nickl on 1/8/19.
//

#include "util.hpp"
#include <algorithm>

namespace MeshSimpl {

bool Internal::sort_and_find_intersection(std::vector<idx>::iterator begin0,
                                          std::vector<idx>::iterator end0,
                                          std::vector<idx>::iterator begin1,
                                          std::vector<idx>::iterator end1) {

    std::sort(begin0, end0);
    std::sort(begin1, end1);
    for (auto i0 = begin0, i1 = begin1; i0 != end0 && i1 != end1;) {
        if (*i0 < *i1)
            ++i0;
        else if (*i0 > *i1)
            ++i1;
        else
            return true;
    }

    return false;
}

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
