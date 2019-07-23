//
// Created by nickl on 1/8/19.
//

#include "util.hpp"

namespace MeshSimpl {

vec3d &operator/=(vec3d &lfs, double rhs) {
  for (auto &x : lfs) x /= rhs;
  return lfs;
}

vec3d operator-(const vec3d &a, const vec3d &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

}  // namespace MeshSimpl
