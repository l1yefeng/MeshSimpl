//
// Created by nickl on 1/8/19.
//

#include "util.hpp"
#include <algorithm>

namespace MeshSimpl {

// Reference: en.cppreference.com/w/cpp/algorithm/set_intersection.html
bool Internal::sort_and_find_intersection(std::vector<idx>& values0,
                                          std::vector<idx>& values1) {
  std::sort(values0.begin(), values0.end());
  std::sort(values1.begin(), values1.end());
  for (auto i0 = values0.begin(), i1 = values1.begin();
       i0 != values0.end() && i1 != values1.end();) {
    if (*i0 < *i1)
      ++i0;
    else if (*i0 > *i1)
      ++i1;
    else
      return true;
  }

  return false;
}

Internal::Quadric& operator+=(Internal::Quadric& lfs,
                              const Internal::Quadric& rhs) {
  for (unsigned int i = 0; i < 10; ++i) lfs[i] += rhs[i];
  return lfs;
}

Internal::Quadric& operator*=(Internal::Quadric& lfs, double rhs) {
  for (auto& x : lfs) x *= rhs;
  return lfs;
}

Internal::Quadric operator+(const Internal::Quadric& a,
                            const Internal::Quadric& b) {
  Internal::Quadric q;
  for (unsigned int i = 0; i < 10; ++i) q[i] = a[i] + b[i];
  return q;
}

vec3d& operator/=(vec3d& lfs, double rhs) {
  for (auto& x : lfs) x /= rhs;
  return lfs;
}

vec3d operator-(const vec3d& a, const vec3d& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

}  // namespace MeshSimpl
