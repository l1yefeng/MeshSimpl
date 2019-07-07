//
// Created by nickl on 5/26/19.
//

#include <initializer_list>

#include "edge.hpp"
#include "faces.hpp"
#include "util.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

bool Faces::onBoundary(idx f) const {
  assert(exists(f));
  for (order k : {0, 1, 2})
    if (side(f, k)->onBoundary()) return true;
  return false;
}

vec3d Faces::vPos(idx f, order k, const Vertices& vertices) const {
  assert(exists(f));
  return vertices[v(f, k)];
}

vec3d Faces::edgeVec(idx f, order k, const Vertices& vertices) const {
  assert(exists(f));
  return vertices[v(f, next(k))] - vertices[v(f, prev(k))];
}

void Faces::compactIndicesAndDie(Indices& indices) {
  for (int lo = 0, hi = size() - 1; true; ++lo, --hi) {
    while (lo <= hi && exists(lo)) ++lo;
    while (lo < hi && !exists(hi)) --hi;
    if (lo >= hi) {
      _indices.resize(lo);
      break;
    }
    std::swap(_indices[lo], _indices[hi]);
  }

  // move to output
  indices = std::move(_indices);
}

}  // namespace Internal
}  // namespace MeshSimpl
