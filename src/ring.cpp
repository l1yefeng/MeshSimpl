//
// Created by nickl on 5/13/19.
//

#include "ring.hpp"
#include <algorithm>     // for sort
#include <cassert>       // for assert
#include "ecol.hpp"      // for update_error_and_center, is_face_elongated
#include "qem_heap.hpp"  // for QEMHeap
#include "vertices.hpp"  // for Vertices

namespace MeshSimpl {
namespace Internal {

bool Ring::check_topo() {
  std::vector<const Neighbor *> star_del, star_kept;
  star_del.reserve(v_del_neighbors.size());
  star_kept.reserve(v_kept_neighbors.size());

  for (const auto &nb : v_del_neighbors) star_del.push_back(&nb);
  for (const auto &nb : v_kept_neighbors) star_kept.push_back(&nb);

  const auto cmp = [&](const Neighbor *nb0, const Neighbor *nb1) -> bool {
    return nb0->second_v(faces) < nb1->second_v(faces);
  };
  std::sort(star_del.begin(), star_del.end(), cmp);
  std::sort(star_kept.begin(), star_kept.end(), cmp);
  for (auto itd = star_del.begin(), itk = star_kept.begin();
       itd != star_del.end() && itk != star_kept.end();) {
    if (cmp(*itd, *itk)) {
      ++itd;
    } else if (cmp(*itk, *itd)) {
      ++itk;
    } else {
      return false;
    }
  }

  return true;
}

bool Ring::check_geom(double foldover_angle) const {
  for (const auto &nb : v_del_neighbors) {
    if (is_face_folded(vertices, faces, nb.f(), nb.center(), edge.col_center(),
                       foldover_angle))
      return false;
  }
  for (const auto &nb : v_kept_neighbors) {
    if (is_face_folded(vertices, faces, nb.f(), nb.center(), edge.col_center(),
                       foldover_angle))
      return false;
  }
  return true;
}

bool Ring::check_quality(double aspect_ratio) const {
  for (const auto &nb : v_del_neighbors) {
    if (is_face_elongated(vertices[v_del], vertices[nb.first_v(faces)],
                          vertices[nb.second_v(faces)], aspect_ratio))
      return false;
  }
  for (const auto &nb : v_del_neighbors) {
    if (is_face_elongated(vertices[v_kept], vertices[nb.first_v(faces)],
                          vertices[nb.second_v(faces)], aspect_ratio))
      return false;
  }
  return true;
}

void InteriorRing::collect() {
  const idx f0 = edge.face(0);
  const idx f1 = edge.face(1);

  Neighbor nb(f0, edge.ord_in_face(0), ccw);
  while (true) {
    nb.rotate(faces);
    if (nb.f() == f1) break;
    v_del_neighbors.push_back(nb);
  }

  bool boundary_hit = false;  // useful in case v_kept is on boundary
  nb = Neighbor(f1, edge.ord_in_face(1), ccw);
  while (true) {
    if (nb.second_edge(faces)->on_boundary()) {
      if (boundary_hit)
        break;  // while-loop breaks here if one v_kept is on boundary
      boundary_hit = true;
      nb = Neighbor(f0, edge.ord_in_face(0), !ccw);
      continue;
    }

    nb.rotate(faces);

    // this face is deleted thus unnecessary to check fold-over if f == f0
    if (nb.f() == f0) {
      // while-loop breaks here if there is no business of boundary
      assert(!boundary_hit);
      break;
    }

    v_kept_neighbors.push_back(nb);
  }
}

bool InteriorRing::check_env() {
  // special case: there are 2 faces, 3 vertices in current component
  // this happens if input contains such component because
  // this edge collapse implementation avoid generating them
  if (v_del_neighbors.empty()) {
    assert(faces[edge.face(0)][edge.ord_in_face(0)] ==
           faces[edge.face(1)][edge.ord_in_face(1)]);
    return false;
  }

  // special case: avoid collapsing a tetrahedron into folded faces
  if (v_del_neighbors.size() == 1 && v_kept_neighbors.size() == 1) return false;

  // fine target to collapse
  return true;
}

void InteriorRing::collapse(QEMHeap &heap, bool fix_boundary) {
  const idx f0 = edge.face(0);
  const idx f1 = edge.face(1);

  // update vertex position and quadric
  if (edge.neither_v_on_border())
    vertices.setPosition(v_kept, edge.col_center());
  vertices.setQ(v_kept, edge.col_q());

  // two kept edges in two deleted faces
  Edge *edge_kept0 = faces.edgeAcrossFrom(f0, v_del);
  Edge *edge_kept1 = faces.edgeAcrossFrom(f1, v_del);

  auto it = v_del_neighbors.begin();
  Edge *dirty_edge_ptr = edge_kept0;

  // first face to process: deleted face 0
  faces.setV(it->f(), it->center(), v_kept);
  heap.erase(faces.side(it->f(), it->j()));

  faces.setSide(it->f(), it->j(), edge_kept0);
  dirty_edge_ptr->replace_f(f0, it->f(), it->j());
  faces.erase(f0);

  // every face centered around deleted vertex
  for (++it; it != v_del_neighbors.end(); ++it) {
    faces.setV(it->f(), it->center(), v_kept);
    dirty_edge_ptr = faces.side(it->f(), it->j());

    dirty_edge_ptr->replace_v(v_del, v_kept, edge.one_v_on_border());

    update_error_and_center(vertices, heap, dirty_edge_ptr, fix_boundary);
  }

  // the deleted face 1
  --it;
  heap.erase(faces.side(it->f(), it->i()));
  faces.setSide(it->f(), it->i(), edge_kept1);
  dirty_edge_ptr = edge_kept1;
  dirty_edge_ptr->replace_f(f1, it->f(), it->i());
  faces.erase(f1);

  vertices.erase(v_del);

  // every edge centered around the kept vertex
  for (auto nb : v_kept_neighbors) {
    dirty_edge_ptr = nb.first_edge(faces);
    if (fix_boundary && dirty_edge_ptr->both_v_on_border()) continue;
    update_error_and_center(vertices, heap, dirty_edge_ptr, fix_boundary);
  }

  if (!(fix_boundary && edge_kept0->both_v_on_border())) {
    update_error_and_center(vertices, heap, edge_kept0, fix_boundary);
  }
}

void BoundaryRing::collect() {
  const idx f = edge.face(0);

  Neighbor nb(f, edge.ord_in_face(0), ccw);
  while (true) {
    const Edge *next_edge = nb.second_edge(faces);

    if (next_edge->on_boundary()) break;

    nb.rotate(faces);
    v_del_neighbors.push_back(nb);
  }

  nb = Neighbor(f, edge.ord_in_face(0), !ccw);
  while (true) {
    const Edge *next_edge = nb.second_edge(faces);

    if (next_edge->on_boundary()) break;

    nb.rotate(faces);
    v_kept_neighbors.push_back(nb);
  }
}

bool BoundaryRing::check_env() {
  // special case: there is ONE triangle in current component
  if (v_kept_neighbors.empty()) return false;

  return true;
}

void BoundaryRing::collapse(QEMHeap &heap, bool fix_boundary) {
  const idx f = edge.face(0);

  // update vertex position and quadric
  vertices.setPosition(v_kept, edge.col_center());
  vertices.setQ(v_kept, edge.col_q());

  Edge *edge_kept = faces.edgeAcrossFrom(f, v_del);
  Edge *dirty_edge_ptr = edge_kept;

  // it is possible that v_del_neighbors is empty
  heap.erase(faces.edgeAcrossFrom(f, v_kept));
  if (v_del_neighbors.empty()) {
    dirty_edge_ptr->drop_f(f);
  } else {
    auto it = v_del_neighbors.begin();
    faces.setSide(it->f(), it->j(), edge_kept);
    dirty_edge_ptr->replace_f(f, it->f(), it->j());
  }

  faces.erase(f);

  // every face centered around deleted vertex
  for (const auto &nb : v_del_neighbors) {
    faces.setV(nb.f(), nb.center(), v_kept);
    dirty_edge_ptr = nb.second_edge(faces);

    dirty_edge_ptr->replace_v(v_del, v_kept, true);

    assert(!dirty_edge_ptr->neither_v_on_border());
    update_error_and_center(vertices, heap, dirty_edge_ptr, false);
  }

  vertices.erase(v_del);

  update_error_and_center(vertices, heap, edge_kept, false);
  for (auto nb : v_kept_neighbors) {
    dirty_edge_ptr = nb.second_edge(faces);
    update_error_and_center(vertices, heap, dirty_edge_ptr, false);
  }
}

}  // namespace Internal
}  // namespace MeshSimpl
