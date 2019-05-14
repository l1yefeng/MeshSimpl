//
// Created by nickl on 5/10/19.
//

#ifndef MESH_SIMPL_MARKER_HPP
#define MESH_SIMPL_MARKER_HPP

#include <cassert>
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Marker {
 private:
  // true means marked deleted for both
  std::vector<bool> vertices;
  std::vector<bool> faces;

 public:
  Marker(size_t vertex_cnt, size_t face_cnt)
      : vertices(vertex_cnt, true), faces(face_cnt, false) {}

  bool exist_f(idx f) const { return !faces[f]; }

  bool exist_v(idx v) const { return !vertices[v]; }

  void mark_v(idx v) {
    assert(!vertices[v]);
    vertices[v] = true;
  }

  void mark_f(idx f) {
    assert(!faces[f]);
    faces[f] = true;
  }

  void mark_unref_v(const std::vector<vec3i> &indices) {
    for (idx f = 0; f < indices.size(); ++f) {
      if (!faces[f]) {
        const auto &vvv = indices[f];
        for (idx v : vvv) vertices[v] = false;
      }
    }
  }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_MARKER_HPP
