//
// Created by nickl on 5/11/19.
//

#include "edge.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

void Edge::attach_1st_face(idx face, order edge_order_in_face) {
  ff[0] = face;
  ord_in_ff[0] = edge_order_in_face;
}

bool Edge::attach_2nd_face(idx face, order edge_order_in_face) {
  assert(!one_v_on_border());
  if (neither_v_on_border()) return false;

  v_on_border[0] = false;
  v_on_border[1] = false;
  ff[1] = face;
  ord_in_ff[1] = edge_order_in_face;
  return true;
}

void Edge::set_v_on_border(bool v0_on_border, bool v1_on_border) {
  v_on_border[0] = v0_on_border;
  v_on_border[1] = v1_on_border;
}

void Edge::plan_collapse(const V &vertices, const Q &quadrics,
                         bool fix_boundary) {
  q = quadrics[vv[0]] + quadrics[vv[1]];

  if (fix_boundary && both_v_on_border()) {
    // the plan is: no plan is needed because it will never by modified
    return;
  } else if (fix_boundary && one_v_on_border()) {
    // the plan is: new position is the position of the vertex on border
    center = vertices[v_on_border[0] ? vv[0] : vv[1]];
    error = q_error(q, center);
  } else {
    // the plan is: new position leads to the lowest error
    const vec3d b{q[6], q[7], q[8]};
    const double c = q[9];

    // computes the inverse of matrix A in quadric
    const double a_det = q[0] * (q[3] * q[5] - q[4] * q[4]) -
                         q[1] * (q[1] * q[5] - q[4] * q[2]) +
                         q[2] * (q[1] * q[4] - q[3] * q[2]);

    if (a_det != 0) {
      // invertible, find position yielding minimal error
      const double a_det_inv = 1.0 / a_det;
      const std::array<double, 6> a_inv{
          (q[3] * q[5] - q[4] * q[4]) * a_det_inv,
          (q[2] * q[4] - q[1] * q[5]) * a_det_inv,
          (q[1] * q[4] - q[2] * q[3]) * a_det_inv,
          (q[0] * q[5] - q[2] * q[2]) * a_det_inv,
          (q[1] * q[2] - q[0] * q[4]) * a_det_inv,
          (q[0] * q[3] - q[1] * q[1]) * a_det_inv,
      };
      center = {-dot({a_inv[0], a_inv[1], a_inv[2]}, b),
                -dot({a_inv[1], a_inv[3], a_inv[4]}, b),
                -dot({a_inv[2], a_inv[4], a_inv[5]}, b)};
      error = dot(b, center) + c;
    } else {
      // not invertible, choose from endpoints and midpoint
      center = midpoint(vertices[vv[0]], vertices[vv[1]]);
      error = q_error(q, center);
      for (const idx v : vv) {
        const double err = q_error(q, vertices[v]);
        if (err < error) {
          center = vertices[v];
          error = err;
        }
      }
    }
  }
}

void Edge::replace_v(idx prev_v, idx new_v, bool new_v_on_border) {
  order ord = v_order(prev_v);
  vv[ord] = new_v;
  if (new_v_on_border) v_on_border[ord] = true;
}

void Edge::drop_f(idx face) {
  assert(!on_boundary());
  order ord = f_order(face);
  ff[ord] = 0;
  ord_in_ff[ord] = INVALID;
  if (ord == 0) swap_faces();

  assert(ord_in_ff[0] != INVALID);
}

void Edge::replace_f(idx prev_f, idx new_f, order new_ord_in_face) {
  order ord = f_order(prev_f);
  assert(ord_in_ff[ord] != INVALID);

  ff[ord] = new_f;
  ord_in_ff[ord] = new_ord_in_face;
}

}  // namespace Internal
}  // namespace MeshSimpl