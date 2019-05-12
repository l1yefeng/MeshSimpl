//
// Created by nickl on 1/8/19.
//

#include "ecol.hpp"

namespace MeshSimpl {
namespace Internal {

static const size_t ESTIMATE_VALENCE = 8;

void update_error_and_center(const V& vertices, const Q& quadrics,
                             QEMHeap& heap, Edge* const edge_ptr,
                             bool fix_boundary) {
  if (fix_boundary && edge_ptr->both_v_on_border()) {
    // this handles the case when non-boundary edge becomes boundary edge cannot
    // erase if given edge was on boundary, since then it would not be in heap
    heap.erase_by_ptr(edge_ptr);
  } else {
    const double error_prev = edge_ptr->col_error();
    edge_ptr->plan_collapse(vertices, quadrics, fix_boundary);
    heap.fix(edge_ptr, error_prev);
  }
}

bool face_fold_over(const V& vertices, const F& indices, const Neighbor& nb,
                    const idx v_move, const vec3d& move_to, double angle) {
  assert(v_move == indices[nb.f()][nb.center()]);
  const idx vi = indices[nb.f()][nb.i()];
  const idx vj = indices[nb.f()][nb.j()];
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

bool extremely_elongated(const V& vertices, const F& indices,
                         const Neighbor& nb, const vec3d& center_pos,
                         double ratio) {
  const vec3d edge_ij =
      vertices[indices[nb.f()][nb.j()]] - vertices[indices[nb.f()][nb.i()]];
  const vec3d edge_ik = center_pos - vertices[indices[nb.f()][nb.i()]];
  const vec3d edge_jk = center_pos - vertices[indices[nb.f()][nb.j()]];
  const double a = magnitude(edge_ij);
  const double b = magnitude(edge_ik);
  const double c = magnitude(edge_jk);
  const double s = (a + b + c) / 2.0;
  const double aspect_ratio = 8.0 * (s - a) * (s - b) * (s - c) / (a * b * c);
  return aspect_ratio < ratio;
}

bool scan_neighbors(const V& vertices, const Connectivity& conn,
                    const Edge& edge, std::vector<Neighbor>& v_del_neighbors,
                    std::vector<idx>& v_kept_neighbor_edges, order del_ord,
                    const SimplifyOptions& options) {
  const idx v_del = edge[del_ord];
  const idx v_kept = edge[1 - del_ord];
  std::vector<Neighbor> v_kept_neighbors;
  std::vector<idx> v_del_twins, v_kept_twins;

  const bool ccw =
      edge.ord_in_face(0) != (conn.v_ord_in_face(edge.face(0), v_del) + 1) % 3;

  if (!edge.on_boundary()) {
    // collapse information around a non-boundary edge

    const idx f0 = edge.face(0);
    const idx f1 = edge.face(1);

    v_del_neighbors.reserve(ESTIMATE_VALENCE);
    v_kept_neighbors.reserve(ESTIMATE_VALENCE);
    v_del_twins.reserve(ESTIMATE_VALENCE);
    v_kept_twins.reserve(ESTIMATE_VALENCE);

    Neighbor nb(f0, edge.ord_in_face(0), ccw);
    while (true) {
      nb.rotate(conn);
      if (nb.f() == f1) break;
      v_del_neighbors.emplace_back(nb);
      v_del_twins.emplace_back(conn.indices[nb.f()][nb.i()]);
    }

    bool boundary_hit = false;  // useful in case v_kept is on boundary
    nb = Neighbor(f1, edge.ord_in_face(1), ccw);
    while (true) {
      const auto& next_edge = conn.edge_of_face(nb.f(), nb.i());

      if (next_edge.on_boundary()) {
        v_kept_twins.emplace_back(conn.indices[nb.f()][nb.j()]);
        v_kept_neighbor_edges.emplace_back(conn.face2edge[nb.f()][nb.i()]);
        if (boundary_hit)
          break;  // while-loop breaks here if one v_kept is on boundary
        boundary_hit = true;
        nb = Neighbor(f0, edge.ord_in_face(0), !ccw);
        continue;
      }

      nb.rotate(conn);
      v_kept_neighbor_edges.emplace_back(conn.face2edge[nb.f()][nb.j()]);
      v_kept_twins.emplace_back(conn.indices[nb.f()][nb.i()]);
      if (nb.f() == f0) {
        // while-loop breaks here if there is no business of boundary
        assert(!boundary_hit);
        break;
      }
      // this face is deleted thus unnecessary to check fold-over if f == f0
      v_kept_neighbors.emplace_back(nb);
    }

    // special case: there are 2 faces, 3 vertices in current component
    // this happens if input contains such component because
    // this edge collapse implementation avoid generating them
    if (v_del_twins.empty()) {
      assert(conn.indices[edge.face(0)][edge.ord_in_face(0)] ==
             conn.indices[edge.face(1)][edge.ord_in_face(1)]);
      return false;
    }

    // we actually want to check intersection
    // on items starting from `v_del_twins[1]` with those in `v_kept_twins`
    v_del_twins.erase(v_del_twins.begin());

    // special case: avoid collapsing a tetrahedron into folded faces
    if (v_del_neighbors.size() == 1 && v_kept_neighbors.size() == 1)
      return false;

  } else {
    // collapse information around a boundary edge

    const idx f = edge.face(0);

    v_del_neighbors.reserve(ESTIMATE_VALENCE / 2);
    v_kept_neighbors.reserve(ESTIMATE_VALENCE / 2);
    v_del_twins.reserve(ESTIMATE_VALENCE / 2);
    v_kept_twins.reserve(ESTIMATE_VALENCE / 2);

    Neighbor nb(f, edge.ord_in_face(0), ccw);
    while (true) {
      const auto& next_edge = conn.edge_of_face(nb.f(), nb.i());

      if (next_edge.on_boundary()) break;

      nb.rotate(conn);
      v_del_neighbors.emplace_back(nb);
      v_del_twins.emplace_back(conn.indices[nb.f()][nb.j()]);
    }

    nb = Neighbor(f, edge.ord_in_face(0), !ccw);
    while (true) {
      const auto& next_edge = conn.edge_of_face(nb.f(), nb.i());
      v_kept_neighbor_edges.emplace_back(conn.face2edge[nb.f()][nb.i()]);

      if (next_edge.on_boundary()) break;

      nb.rotate(conn);
      v_kept_neighbors.emplace_back(nb);
      v_kept_twins.emplace_back(conn.indices[nb.f()][nb.j()]);
    }

    // special case: there is ONE triangle in current component
    if (v_kept_neighbors.empty()) return false;
    if (v_kept_neighbors.empty()) return false;
  }

  // check connectivity: intersection nonempty indicates topology change after
  // collapse
  if (sort_and_find_intersection(v_del_twins, v_kept_twins)) return false;

  // check geometry
  for (const auto& nb_del : v_del_neighbors)
    if (face_fold_over(vertices, conn.indices, nb_del, v_del, edge.col_center(),
                       options.fold_over_angle_threshold))
      return false;
  for (const auto& nb_kept : v_kept_neighbors)
    if (face_fold_over(vertices, conn.indices, nb_kept, v_kept,
                       edge.col_center(), options.fold_over_angle_threshold))
      return false;

  // check face quality
  for (const auto& nb_del : v_del_neighbors)
    if (extremely_elongated(vertices, conn.indices, nb_del, edge.col_center(),
                            options.aspect_ratio_at_least))
      return false;
  for (const auto& nb_kept : v_del_neighbors)
    if (extremely_elongated(vertices, conn.indices, nb_kept, edge.col_center(),
                            options.aspect_ratio_at_least))
      return false;

  // good to go
  return true;
}

bool edge_collapse(V& vertices, Internal::Connectivity& conn, Q& quadrics,
                   QEMHeap& heap, const idx ecol_target, order del_ord,
                   const SimplifyOptions& options) {
  auto& edge = conn.edges[ecol_target];

  // if non-boundary edge has two endpoints on boundary, we avoid collapsing it
  // because it is possible to produce non-manifold vertex
  if (edge.both_v_on_border() && !edge.on_boundary()) {
    heap.penalize(ecol_target);
    return false;
  }

  // collect neighboring faces into two containers
  std::vector<Neighbor> v_del_neighbors;
  std::vector<idx> v_kept_neighbor_edges;
  if (!scan_neighbors(vertices, conn, edge, v_del_neighbors,
                      v_kept_neighbor_edges, del_ord, options)) {
    heap.penalize(ecol_target);
    return false;
  }

  // now that we are certain this edge is to be collapsed, remove it from heap
  heap.pop();

  const idx v_del = edge[del_ord];
  const idx v_kept = edge[1 - del_ord];

  if (!edge.on_boundary()) {
    // collapse a non-boundary edge and update its neighborhood

    const idx f0 = edge.face(0);
    const idx f1 = edge.face(1);

    // update vertex position and quadric
    if (edge.neither_v_on_border()) vertices[v_kept] = edge.col_center();
    quadrics[v_kept] = edge.col_q();

    vec2i e_kept;  // global index of kept edge in two deleted faces
    e_kept[0] = conn.edge_idx_across_from_v(f0, v_del);
    e_kept[1] = conn.edge_idx_across_from_v(f1, v_del);

    auto it = v_del_neighbors.begin();
    Edge* dirty_edge_ptr = &conn.edges[e_kept[0]];

    // first face to process: deleted face 0
    conn.indices[it->f()][it->center()] = v_kept;
    heap.erase(conn.face2edge[it->f()][it->j()]);
    conn.face2edge[it->f()][it->j()] = e_kept[0];
    dirty_edge_ptr->replace_f(f0, it->f(), it->j());

    // every face centered around deleted vertex
    for (++it; it != v_del_neighbors.end(); ++it) {
      conn.indices[it->f()][it->center()] = v_kept;
      dirty_edge_ptr = &conn.edge_of_face(it->f(), it->j());

      dirty_edge_ptr->replace_v(v_del, v_kept, edge.one_v_on_border());

      update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr,
                              options.fix_boundary);
    }

    // the deleted face 1
    --it;
    heap.erase(conn.face2edge[it->f()][it->i()]);
    conn.face2edge[it->f()][it->i()] = e_kept[1];
    dirty_edge_ptr = &conn.edges[e_kept[1]];
    dirty_edge_ptr->replace_f(f1, it->f(), it->i());

    // every edge centered around the kept vertex
    for (const idx e_dirty : v_kept_neighbor_edges)
      if (!(options.fix_boundary && conn.edges[e_dirty].both_v_on_border()))
        update_error_and_center(vertices, quadrics, heap, &conn.edges[e_dirty],
                                options.fix_boundary);

  } else {
    // collapse a boundary edge and update its neighborhood

    const idx f = edge.face(0);

    // update vertex position and quadric
    vertices[v_kept] = edge.col_center();
    quadrics[v_kept] = edge.col_q();

    const idx e_kept_idx = conn.edge_idx_across_from_v(f, v_del);
    Edge* dirty_edge_ptr = &conn.edges[e_kept_idx];

    // it is possible that v_del_neighbors is empty
    heap.erase(conn.edge_idx_across_from_v(f, v_kept));
    if (v_del_neighbors.empty()) {
      dirty_edge_ptr->drop_f(f);
    } else {
      auto it = v_del_neighbors.begin();
      conn.face2edge[it->f()][it->j()] = e_kept_idx;
      dirty_edge_ptr->replace_f(f, it->f(), it->j());
    }

    // every face centered around deleted vertex
    for (const auto& nb : v_del_neighbors) {
      conn.indices[nb.f()][nb.center()] = v_kept;
      dirty_edge_ptr = &conn.edge_of_face(nb.f(), nb.i());

      dirty_edge_ptr->replace_v(v_del, v_kept, true);

      assert(!dirty_edge_ptr->neither_v_on_border());
      update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr, false);
    }

    assert(!v_kept_neighbor_edges.empty());
    for (const idx e_dirty : v_kept_neighbor_edges)
      update_error_and_center(vertices, quadrics, heap, &conn.edges[e_dirty],
                              false);
  }

  return true;
}

}  // namespace Internal
}  // namespace MeshSimpl
