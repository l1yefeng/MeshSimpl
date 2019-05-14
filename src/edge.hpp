//
// Created by nickl on 03/20/19.
//

#ifndef MESH_SIMPL_EDGE_HPP
#define MESH_SIMPL_EDGE_HPP

#include <cassert>
#include <limits>
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

// Edge defines the struct of an edge
//
// Only some convenient methods are defined for the purpose of less repetition.
// Methods to initialize (adding first and second face), and update are not
// provided.
class Edge {
 private:
  static const order INVALID = -1;

  // sum of quadrics of two endpoints
  Quadric q = {};

  // where this edge collapse into
  vec3d center = {};

  // quadric error value
  double error = 0;

  // index of two endpoints
  vec2i vv;

  // index of two incident faces;
  // boundary edges have faces[0] and ord_in_face[1] is INVALID;
  // all edges have at least face[0] after edges are constructed in
  // `construct_edges()`
  vec2i ff = {};

  // order of this edge in two incident faces, match the order to member "ff"
  std::array<order, 2> ord_in_ff = {INVALID, INVALID};

  // vertices on boundary
  std::array<bool, 2> v_on_border;

 public:
  // Construct an edge with two endpoints indexes.
  // After constructed this edge
  //  - (Required) call attach_1st_face()
  //  - (Optional) call attach_2nd_face()
  //  - (Required) call set_v_on_border()
  // Then this edge is correctly initialized
  Edge(idx v0, idx v1) : vv({v0, v1}), v_on_border({true, true}) {}

  // Attach a face to edge.
  // No other modifier member functions should be called before this call
  void attach_1st_face(idx face, order edge_order_in_face);

  // Attach a second face to edge.
  // Returns false if this is not the first time attach_2nd_face() being called.
  // Must be called after attach_1st_face();
  // No other modifier member functions should be called before this call
  bool attach_2nd_face(idx face, order edge_order_in_face);

  // Before this call, the vertex-on-border information is
  //  - both_on_border: this edge is on border;
  //  - neither_on_border: this edge is not on border;
  // Only after this call does the vertex-on-border information actually mean
  // endpoint on border or not
  void set_v_on_border(bool v0_on_border, bool v1_on_border);

  //
  // public methods for retrieval of information
  //

  // Returns the collapse center of next collapse
  const vec3d &col_center() const { return center; }

  // Returns the error value made of next collapse
  double col_error() const { return error; }

  // Returns the sum of quadrics of two endpoints
  const Quadric &col_q() const { return q; }

  bool both_v_on_border() const { return v_on_border[0] && v_on_border[1]; }

  bool neither_v_on_border() const {
    return !v_on_border[0] && !v_on_border[1];
  }

  bool one_v_on_border() const {
    return !both_v_on_border() && !neither_v_on_border();
  }

  idx face(order ord) const { return ff[ord]; }

  // Returns the order of this edge in face(ord)
  order ord_in_face(order ord) const { return ord_in_ff[ord]; }

  // Returns true if egde is on border
  // NOTE: different from both endpoints on border
  bool on_boundary() const { return ord_in_ff[1] == INVALID; }

  // Returns the order of face f on this edge
  order f_order(idx f) const {
    assert(ff[0] == f || (ff[1] == f && ord_in_ff[1] != INVALID));
    return ff[0] == f ? 0 : 1;
  }

  // Returns the order of endpoint v on this edge
  order v_order(idx v) const {
    assert(vv[0] == v || vv[1] == v);
    return vv[0] == v ? 0 : 1;
  }

  // Returns the order of the endpoint next collapse should delete.
  // By keeping one endpoint and deleting the other, followed by update on the
  // kept endpoint (position, quadric, etc), the edge is collapsed into a new
  // vertex (the center)
  order v_del_order() const { return v_on_border[0] ? 1 : 0; }

  const vec2i &faces() const { return ff; }

  const vec2i &vertices() const { return vv; }

  idx operator[](order ord) const { return vv[ord]; }

  //
  // public methods for update
  //

  // Plane next collapse.
  // Will set
  //  - what is current sum of endpoints quadrics
  //  - which position to collapse into (center)
  //  - what will be the error
  void plan_collapse(const V &vertices, const Q &quadrics, bool fix_boundary);

  void set_infty_error() { error = std::numeric_limits<double>::max(); }

  void swap_faces() {
    std::swap(ff[0], ff[1]);
    std::swap(ord_in_ff[0], ord_in_ff[1]);
  }

  void replace_v(idx prev_v, idx new_v, bool new_v_on_border);

  void drop_f(idx face);

  void replace_f(idx prev_f, idx new_f, order new_ord_in_face);
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_EDGE_HPP
