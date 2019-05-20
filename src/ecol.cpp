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

bool is_face_folded(const V &vertices, const F &faces, const Face &face,
                    order moved, const vec3d &position, double angle) {
  const idx vk = face[moved];
  const idx vi = face[next(moved)];
  const idx vj = face[prev(moved)];
  const vec3d vec_ij = vertices[vj] - vertices[vi];
  const vec3d vec_jk_prv = vertices[vk] - vertices[vj];
  const vec3d vec_jk_new = position - vertices[vj];
  vec3d normal_prv = cross(vec_ij, vec_jk_prv);
  double mag_prv = magnitude(normal_prv);
  assert(mag_prv != 0);
  vec3d normal_new = cross(vec_ij, vec_jk_new);
  double mag_new = magnitude(normal_new);
  if (mag_new == 0) return true;
  normal_prv /= mag_prv;
  normal_new /= mag_new;
  double cos = dot(normal_prv, normal_new);
  return cos < angle;
}

bool is_face_elongated(const vec3d &pos0, const vec3d &pos1, const vec3d &pos2,
                       double ratio) {
  const vec3d vec01 = pos1 - pos0;
  const vec3d vec02 = pos2 - pos0;
  const vec3d vec12 = pos2 - pos1;
  const double a = magnitude(vec01);
  const double b = magnitude(vec12);
  const double c = magnitude(vec02);
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
