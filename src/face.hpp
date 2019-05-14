//
// Created by nickl on 5/12/19.
//

#ifndef MESH_SIMPL_FACE_HPP
#define MESH_SIMPL_FACE_HPP

#include "edge.hpp"
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Face {
 private:
  vec3i &vvv;
  std::array<Edge *, 3> eee = {};

 public:
  explicit Face(vec3i &vertices) : vvv(vertices) {}

  Edge *edge(order ord) const { return eee[ord]; }

  void attach_edge(Edge *edge, order ord) { eee[ord] = edge; }

  // accessors

  Edge *edge_across_from(idx v) const { return eee[v_order(v)]; }

  order v_order(idx v) const {
    if (v == vvv[0])
      return 0;
    else if (v == vvv[1])
      return 1;
    else if (v == vvv[2])
      return 2;
    else
      assert(false);
    return -1;
  }

  const vec3i &vertices() const { return vvv; }

  idx operator[](order ord) const { return vvv[ord]; }

  // modifiers

  void replace_v(order ord, idx new_v) { vvv[ord] = new_v; }

  void replace_edge(order ord, Edge *new_edge) { eee[ord] = new_edge; }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_FACE_HPP
