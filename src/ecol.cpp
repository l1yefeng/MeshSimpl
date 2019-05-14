//
// Created by nickl on 1/8/19.
//

#include "ecol.hpp"
#include <memory>
#include "ring.hpp"

namespace MeshSimpl {
namespace Internal {

void update_error_and_center(const V &vertices, const Q &quadrics,
                             QEMHeap &heap, Edge *const edge_ptr,
                             bool fix_boundary) {
  if (fix_boundary && edge_ptr->both_v_on_border()) {
    // this handles the case when non-boundary edge becomes boundary edge cannot
    // erase if given edge was on boundary, since then it would not be in heap
    heap.erase(edge_ptr);
  } else {
    const double error_prev = edge_ptr->col_error();
    edge_ptr->plan_collapse(vertices, quadrics, fix_boundary);
    heap.fix(edge_ptr, error_prev);
  }
}

bool face_fold_over(const V &vertices, const F &faces, const Neighbor &nb,
                    idx v_move, const vec3d &move_to, double angle) {
  assert(v_move == faces[nb.f()][nb.center()]);
  const idx vi = faces[nb.f()][nb.i()];
  const idx vj = faces[nb.f()][nb.j()];
  const vec3d e0 = vertices[vj] - vertices[vi];
  const vec3d e1_prev = vertices[v_move] - vertices[vj];
  const vec3d e1_new = move_to - vertices[vj];
  vec3d normal_prev = cross(e0, e1_prev);
  double normal_prev_mag = magnitude(normal_prev);
  vec3d normal_new = cross(e0, e1_new);
  double normal_new_mag = magnitude(normal_new);
  if (normal_new_mag == 0) return true;
  normal_prev /= normal_prev_mag;
  normal_new /= normal_new_mag;
  double cos = dot(normal_prev, normal_new);
  return cos < angle;
}

bool extremely_elongated(const V &vertices, const F &faces, const Neighbor &nb,
                         const vec3d &center_pos, double ratio) {
  const vec3d edge_ij =
      vertices[faces[nb.f()][nb.j()]] - vertices[faces[nb.f()][nb.i()]];
  const vec3d edge_ik = center_pos - vertices[faces[nb.f()][nb.i()]];
  const vec3d edge_jk = center_pos - vertices[faces[nb.f()][nb.j()]];
  const double a = magnitude(edge_ij);
  const double b = magnitude(edge_ik);
  const double c = magnitude(edge_jk);
  const double s = (a + b + c) / 2.0;
  const double aspect_ratio = 8.0 * (s - a) * (s - b) * (s - c) / (a * b * c);
  return aspect_ratio < ratio;
}

bool edge_collapse(V &vertices, F &faces, Q &quadrics, QEMHeap &heap,
                   Edge &target, const SimplifyOptions &options) {
  // if non-boundary edge has two endpoints on boundary, we avoid collapsing it
  // because it is possible to produce non-manifold vertex
  if (target.both_v_on_border() && !target.on_boundary()) {
    heap.penalize(&target);
    return false;
  }

  std::unique_ptr<Ring> ring;
  if (target.on_boundary())
    ring = std::unique_ptr<Ring>(new BoundaryRing(vertices, faces, target));
  else
    ring = std::unique_ptr<Ring>(new InteriorRing(vertices, faces, target));

  ring->collect();

  if (!(ring->check_env() && ring->check_topo() &&
        ring->check_geom(options.fold_over_angle_threshold) &&
        ring->check_quality(options.aspect_ratio_at_least))) {
    heap.penalize(&target);
    return false;
  }

  // now that we are certain this edge is to be collapsed, remove it from heap
  heap.pop();

  ring->collapse(quadrics, heap, options.fix_boundary);

  return true;
}

}  // namespace Internal
}  // namespace MeshSimpl
