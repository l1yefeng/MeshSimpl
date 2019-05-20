//
// Created by nickl on 3/13/19.
//

#ifndef MESH_SIMPL_NEIGHBOR_HPP
#define MESH_SIMPL_NEIGHBOR_HPP

#include "edge.hpp"
#include "face.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

// This is a compact structure to represent an incident face around a center
// vertex. It is created for the traverse of the incident faces of a vertex
// (center) during edge-collapse operation.
//
// It contains a face index and two vertex (the third is the center) order: i
// and j. i, j, and center are local orders of three vertices. In typical
// (clockwise) cases, the orientation is center -> i -> j -> center
//
// Mind that in a face, any edge and the vertex across from it have the same
// order
class Neighbor {
 private:
  idx face;  // face index in indices
  bool ccw;  // determines the direction of rotation
  order vi;  // i = (j + 1) % 3 if clockwise
  order vj;  // j = (center + 1) % 3 if clockwise

  order get_i_from_j(order j) { return ccw ? prev(j) : next(j); }

 public:
  Neighbor(idx face, order vj, bool ccw)
      : face(face), ccw(ccw), vi(get_i_from_j(vj)), vj(vj) {}

  idx f() const { return face; }

  order i() const { return vi; }

  order j() const { return vj; }

  order center() const { return 3 - vi - vj; }

  void rotate(const F& faces) {
    const Edge* curr_edge = second_edge(faces);
    assert(!curr_edge->on_boundary());

    const idx prev_center = faces[face][center()];

    const order ford = curr_edge->f_order(face);
    const idx next_face = curr_edge->face(1 - ford);
    assert(face != next_face);
    face = next_face;
    vj = curr_edge->ord_in_face(1 - ford);
    vi = get_i_from_j(vj);

    assert(faces[face][center()] == prev_center);
  }

  Edge* first_edge(const F& faces) const { return faces[f()].edge(vj); }

  Edge* second_edge(const F& faces) const { return faces[f()].edge(vi); }

  idx first_v(const F& faces) const { return faces[f()][vi]; }

  idx second_v(const F& faces) const { return faces[f()][vj]; }

  const Face& f(const F& faces) const { return faces[face]; }

  idx center_v(const F& faces) const { return f(faces)[center()]; }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_NEIGHBOR_HPP
