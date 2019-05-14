//
// Created by nickl on 5/13/19.
//

#include "ring.hpp"
#include "ecol.hpp"

namespace MeshSimpl {
namespace Internal {

bool Ring::check_geom(double foldover_angle) const {
  for (const auto &nb_del : v_del_neighbors) {
    if (face_fold_over(vertices, faces, nb_del, v_del, edge.col_center(),
                       foldover_angle))
      return false;
  }
  for (const auto &nb_kept : v_kept_neighbors) {
    if (face_fold_over(vertices, faces, nb_kept, v_kept, edge.col_center(),
                       foldover_angle))
      return false;
  }
  return true;
}

bool Ring::check_quality(double aspect_ratio) const {
  for (const auto &nb_del : v_del_neighbors) {
    if (extremely_elongated(vertices, faces, nb_del, edge.col_center(),
                            aspect_ratio))
      return false;
  }
  for (const auto &nb_kept : v_del_neighbors) {
    if (extremely_elongated(vertices, faces, nb_kept, edge.col_center(),
                            aspect_ratio))
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
    v_del_twins.push_back(faces[nb.f()][nb.i()]);
  }

  bool boundary_hit = false;  // useful in case v_kept is on boundary
  nb = Neighbor(f1, edge.ord_in_face(1), ccw);
  while (true) {
    const Edge *next_edge = faces[nb.f()].edge(nb.i());

    if (next_edge->on_boundary()) {
      v_kept_twins.push_back(faces[nb.f()][nb.j()]);
      v_kept_neighbor_edges.push_back(faces[nb.f()].edge(nb.i()));
      if (boundary_hit)
        break;  // while-loop breaks here if one v_kept is on boundary
      boundary_hit = true;
      nb = Neighbor(f0, edge.ord_in_face(0), !ccw);
      continue;
    }

    nb.rotate(faces);
    v_kept_neighbor_edges.push_back(faces[nb.f()].edge(nb.j()));
    v_kept_twins.push_back(faces[nb.f()][nb.i()]);

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
  if (v_del_twins.empty()) {
    assert(faces[edge.face(0)][edge.ord_in_face(0)] ==
           faces[edge.face(1)][edge.ord_in_face(1)]);
    return false;
  }

  // we actually want to check intersection
  // on items starting from `v_del_twins[1]` with those in `v_kept_twins`
  v_del_twins.erase(v_del_twins.begin());

  // special case: avoid collapsing a tetrahedron into folded faces
  if (v_del_neighbors.size() == 1 && v_kept_neighbors.size() == 1) return false;

  // fine target to collapse
  return true;
}

void InteriorRing::collapse(Q &quadrics, QEMHeap &heap, bool fix_boundary) {
  const idx f0 = edge.face(0);
  const idx f1 = edge.face(1);

  // update vertex position and quadric
  if (edge.neither_v_on_border()) vertices[v_kept] = edge.col_center();
  quadrics[v_kept] = edge.col_q();

  // two kept edges in two deleted faces
  Edge *edge_kept0 = faces[f0].edge_across_from(v_del);
  Edge *edge_kept1 = faces[f1].edge_across_from(v_del);

  auto it = v_del_neighbors.begin();
  Edge *dirty_edge_ptr = edge_kept0;

  // first face to process: deleted face 0
  faces[it->f()].replace_v(it->center(), v_kept);
  heap.erase(faces[it->f()].edge(it->j()));

  faces[it->f()].replace_edge(it->j(), edge_kept0);
  dirty_edge_ptr->replace_f(f0, it->f(), it->j());

  // every face centered around deleted vertex
  for (++it; it != v_del_neighbors.end(); ++it) {
    faces[it->f()].replace_v(it->center(), v_kept);
    dirty_edge_ptr = faces[it->f()].edge(it->j());

    dirty_edge_ptr->replace_v(v_del, v_kept, edge.one_v_on_border());

    update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr,
                            fix_boundary);
  }

  // the deleted face 1
  --it;
  heap.erase(faces[it->f()].edge(it->i()));
  faces[it->f()].replace_edge(it->i(), edge_kept1);
  dirty_edge_ptr = edge_kept1;
  dirty_edge_ptr->replace_f(f1, it->f(), it->i());

  // every edge centered around the kept vertex
  for (Edge *e_dirty : v_kept_neighbor_edges) {
    if (fix_boundary && e_dirty->both_v_on_border()) continue;
    update_error_and_center(vertices, quadrics, heap, e_dirty, fix_boundary);
  }
}

void BoundaryRing::collect() {
  const idx f = edge.face(0);

  Neighbor nb(f, edge.ord_in_face(0), ccw);
  while (true) {
    const Edge *next_edge = faces[nb.f()].edge(nb.i());

    if (next_edge->on_boundary()) break;

    nb.rotate(faces);
    v_del_neighbors.push_back(nb);
    v_del_twins.push_back(faces[nb.f()][nb.j()]);
  }

  nb = Neighbor(f, edge.ord_in_face(0), !ccw);
  while (true) {
    const Edge *next_edge = faces[nb.f()].edge(nb.i());
    v_kept_neighbor_edges.push_back(faces[nb.f()].edge(nb.i()));

    if (next_edge->on_boundary()) break;

    nb.rotate(faces);
    v_kept_neighbors.push_back(nb);
    v_kept_twins.push_back(faces[nb.f()][nb.j()]);
  }
}

bool BoundaryRing::check_env() {
  // special case: there is ONE triangle in current component
  if (v_kept_neighbors.empty()) return false;

  return true;
}

void BoundaryRing::collapse(Q &quadrics, QEMHeap &heap, bool) {
  const idx f = edge.face(0);

  // update vertex position and quadric
  vertices[v_kept] = edge.col_center();
  quadrics[v_kept] = edge.col_q();

  Edge *edge_kept = faces[f].edge_across_from(v_del);
  Edge *dirty_edge_ptr = edge_kept;

  // it is possible that v_del_neighbors is empty
  heap.erase(faces[f].edge_across_from(v_kept));
  if (v_del_neighbors.empty()) {
    dirty_edge_ptr->drop_f(f);
  } else {
    auto it = v_del_neighbors.begin();
    faces[it->f()].replace_edge(it->j(), edge_kept);
    dirty_edge_ptr->replace_f(f, it->f(), it->j());
  }

  // every face centered around deleted vertex
  for (const auto &nb : v_del_neighbors) {
    faces[nb.f()].replace_v(nb.center(), v_kept);
    dirty_edge_ptr = faces[nb.f()].edge(nb.i());

    dirty_edge_ptr->replace_v(v_del, v_kept, true);

    assert(!dirty_edge_ptr->neither_v_on_border());
    update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr, false);
  }

  assert(!v_kept_neighbor_edges.empty());
  for (Edge *e_dirty : v_kept_neighbor_edges) {
    update_error_and_center(vertices, quadrics, heap, e_dirty, false);
  }
}

}  // namespace Internal
}  // namespace MeshSimpl
