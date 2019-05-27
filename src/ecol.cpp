//
// Created by nickl on 1/8/19.
//

#include "ecol.hpp"
#include <cassert>       // for assert
#include <memory>        // for unique_ptr
#include "edge.hpp"      // for Edge
#include "faces.hpp"     // for Faces
#include "qem_heap.hpp"  // for QEMHeap
#include "ring.hpp"      // for Ring, BoundaryRing, InteriorRing
#include "util.hpp"      // for operator-, magnitude, cross, operator/=, dot

namespace MeshSimpl {
namespace Internal {

void update_error_and_center(Vertices &vertices, QEMHeap &heap,
                             Edge *const edge_ptr, bool fix_boundary) {
  if (fix_boundary && edge_ptr->both_v_on_border(vertices)) {
    // this handles the case when non-boundary edge becomes boundary edge cannot
    // erase if given edge was on boundary, since then it would not be in heap
    heap.erase(edge_ptr);
  } else {
    const double error_prev = edge_ptr->col_error();
    edge_ptr->plan_collapse(vertices, fix_boundary);
    heap.fix(edge_ptr, error_prev);
  }
}

bool is_face_folded(const Vertices &vertices, const Faces &faces, idx f,
                    order moved, const vec3d &position, double angle) {
  const vec3d &vk = faces.vPos(f, moved, vertices);
  const vec3d &vi = faces.vPos(f, next(moved), vertices);
  const vec3d &vj = faces.vPos(f, prev(moved), vertices);
  const vec3d edgeVec0 = vj - vi;
  const vec3d edgeVec1 = vk - vj;
  const vec3d edgeVec1New = position - vj;
  vec3d normal_prv = cross(edgeVec0, edgeVec1);
  double mag_prv = magnitude(normal_prv);
  assert(mag_prv != 0);
  vec3d normal_new = cross(edgeVec0, edgeVec1New);
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

int edge_collapse(Vertices &vertices, Faces &faces, QEMHeap &heap, Edge &target,
                  const SimplifyOptions &options) {
  // if non-boundary edge has two endpoints on boundary, we avoid collapsing it
  // because it is possible to produce non-manifold vertex
  if (target.both_v_on_border(vertices) && !target.on_boundary()) {
    heap.penalize(&target);
    return 0;
  }

  std::unique_ptr<Ring> ring;
  if (target.on_boundary())
    ring = std::unique_ptr<Ring>(new BoundaryRing(vertices, faces, target));
  else
    ring = std::unique_ptr<Ring>(new InteriorRing(vertices, faces, target));

  ring->collect();

  if (!ring->check(options)) {
    heap.penalize(&target);
    return 0;
  }

  // now that we are certain this edge is to be collapsed, remove it from heap
  heap.pop();

  ring->collapse(heap, options.fix_boundary);

  return 1;
}

}  // namespace Internal
}  // namespace MeshSimpl
