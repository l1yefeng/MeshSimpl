//
// Created by nickl on 5/25/19.
//

#ifndef MESH_SIMPL_FACES_HPP
#define MESH_SIMPL_FACES_HPP

#include <array>
#include <cassert>
#include <utility>
#include <vector>
#include <ext/alloc_traits.h>

#include "erasable.hpp"
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Vertices;
class Edge;

class Faces : public Erasables {
 private:
  Indices _indices;
  // FIXME: BAD idea to use pointers. this program never resizes the vector
  // edges so these pointers are never invalidated, but it still is bad manner
  std::vector<std::array<Edge*, 3>> _sides;

 public:
  // Embed indices and allocate space for sides
  explicit Faces(Indices& indices)
      : Erasables(indices.size()),
        _indices(std::move(indices)),
        _sides(size()) {}

  // Get/set side: edge of a face
  Edge* side(idx f, order ord) const {
    assert(exists(f));
    return _sides[f][ord];
  }
  void setSide(idx f, order ord, Edge* edge) {
    assert(exists(f));
    _sides[f][ord] = edge;
  }
  Edge* edgeAcrossFrom(idx f, idx v) const { return side(f, orderOf(f, v)); }

  // Get/set vertex index
  idx v(idx f, order ord) const {
    assert(exists(f));
    return indices(f)[ord];
  }
  void setV(idx f, order ord, idx v) {
    assert(exists(f));
    _indices[f][ord] = v;
  }
  const vec3i& indices(idx f) const {
    assert(exists(f));
    return _indices[f];
  }
  const vec3i& operator[](idx f) const { return indices(f); }

  bool onBoundary(idx f) const;

  vec3d vPos(idx f, order k, const Vertices& vertices) const;

  vec3d edgeVec(idx f, order k, const Vertices& vertices) const;

  void compactIndicesAndDie(Indices& indices);

  // Return the order of a vertex v in some face f
  order orderOf(idx f, idx v) const {
    if (v == indices(f)[0])
      return 0;
    else if (v == indices(f)[1])
      return 1;
    else if (v == indices(f)[2])
      return 2;
    else
      assert(false);
    return INVALID;
  }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_FACES_HPP
