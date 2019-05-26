//
// Created by nickl on 5/26/19.
//

#include "faces.hpp"
#include <initializer_list>  // for initializer_list
#include "edge.hpp"          // for Edge
#include "util.hpp"          // for next, operator-, prev
#include "vertices.hpp"      // for Vertices

namespace MeshSimpl {
namespace Internal {

bool Faces::onBoundary(idx f) const {
  assert(exists(f));
  for (order k : {0, 1, 2})
    if (side(f, k)->on_boundary()) return true;
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

}  // namespace Internal
}  // namespace MeshSimpl
