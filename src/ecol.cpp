//
// Created by nickl on 1/8/19.
//

#include "ecol.hpp"

namespace MeshSimpl {
namespace Internal {

static const size_t ESTIMATE_VALENCE = 8;

void optimal_ecol_vertex_placement(const V& vertices, Edge& edge) {
  const Quadric& q = edge.q;
  const vec3d b{q[6], q[7], q[8]};
  const double c = q[9];

  // computes the inverse of matrix A in quadric
  const double a_det = q[0] * (q[3] * q[5] - q[4] * q[4]) -
                       q[1] * (q[1] * q[5] - q[4] * q[2]) +
                       q[2] * (q[1] * q[4] - q[3] * q[2]);

  if (a_det != 0) {
    // invertible, find position yielding minimal error
    const double a_det_inv = 1.0 / a_det;
    const std::array<double, 6> a_inv{(q[3] * q[5] - q[4] * q[4]) * a_det_inv,
                                      (q[2] * q[4] - q[1] * q[5]) * a_det_inv,
                                      (q[1] * q[4] - q[2] * q[3]) * a_det_inv,
                                      (q[0] * q[5] - q[2] * q[2]) * a_det_inv,
                                      (q[1] * q[2] - q[0] * q[4]) * a_det_inv,
                                      (q[0] * q[3] - q[1] * q[1]) * a_det_inv};
    edge.center = {-dot({a_inv[0], a_inv[1], a_inv[2]}, b),
                   -dot({a_inv[1], a_inv[3], a_inv[4]}, b),
                   -dot({a_inv[2], a_inv[4], a_inv[5]}, b)};
    edge.error = dot(b, edge.center) + c;
  } else {
    // not invertible, choose from endpoints and midpoint
    edge.center =
        midpoint(vertices[edge.vertices[0]], vertices[edge.vertices[1]]);
    edge.error = q_error(edge.q, edge.center);
    for (const idx v : edge.vertices) {
      const double err = q_error(edge.q, vertices[v]);
      if (err < edge.error) {
        edge.center = vertices[v];
        edge.error = err;
      }
    }
  }
}

void set_edge_error(const V& vertices, const Q& quadrics, Edge& edge,
                    bool fix_boundary) {
  // error = v(Q1+Q2)v = vQv, v is new vertex position after edge-collapse
  const auto& vv = edge.vertices;
  edge.q = quadrics[vv[0]] + quadrics[vv[1]];

  if (fix_boundary &&
      (edge.boundary_v == Edge::V0 || edge.boundary_v == Edge::V1)) {
    edge.center = vertices[vv[edge.boundary_v]];
    edge.error = q_error(edge.q, edge.center);
  } else {
    optimal_ecol_vertex_placement(vertices, edge);
  }
}

void compute_errors(const V& vertices, const Q& quadrics, E& edges,
                    bool fix_boundary) {
  for (auto& edge : edges)
    set_edge_error(vertices, quadrics, edge, fix_boundary);
}

void update_error_and_center(const V& vertices, const Q& quadrics,
                             QEMHeap& heap, Edge* const edge_ptr,
                             bool fix_boundary) {
  if (fix_boundary && edge_ptr->boundary_v == Edge::BOTH) {
    // this handles the case when non-boundary edge becomes boundary edge cannot
    // erase if given edge was on boundary, since then it would not be in heap
    heap.erase_by_ptr(edge_ptr);
  } else {
    const double error_prev = edge_ptr->error;
    set_edge_error(vertices, quadrics, *edge_ptr, fix_boundary);
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
                    std::vector<idx>& v_kept_neighbor_edges,
                    const SimplifyOptions& options) {
  const idx v_del = edge.vertices[edge.v_del_order()];
  const idx v_kept = edge.vertices[1 - edge.v_del_order()];
  std::vector<Neighbor> v_kept_neighbors;
  std::vector<idx> v_del_twins, v_kept_twins;

  const bool ccw = edge.ord_in_faces[0] !=
                   (conn.v_ord_in_face(edge.faces[0], v_del) + 1) % 3;

  if (!edge.on_boundary()) {
    // collapse information around a non-boundary edge

    const vec2i& ff = edge.faces;

    v_del_neighbors.reserve(ESTIMATE_VALENCE);
    v_kept_neighbors.reserve(ESTIMATE_VALENCE);
    v_del_twins.reserve(ESTIMATE_VALENCE);
    v_kept_twins.reserve(ESTIMATE_VALENCE);

    Neighbor nb(ff[0], edge.ord_in_faces[0], ccw);
    while (true) {
      nb.rotate(conn);
      if (nb.f() == ff[1]) break;
      v_del_neighbors.emplace_back(nb);
      v_del_twins.emplace_back(conn.indices[nb.f()][nb.i()]);
    }

    bool boundary_hit = false;  // useful in case v_kept is on boundary
    nb = Neighbor(ff[1], edge.ord_in_faces[1], ccw);
    while (true) {
      const auto& next_edge = conn.edge_of_face(nb.f(), nb.i());

      if (next_edge.on_boundary()) {
        v_kept_twins.emplace_back(conn.indices[nb.f()][nb.j()]);
        v_kept_neighbor_edges.emplace_back(conn.face2edge[nb.f()][nb.i()]);
        if (boundary_hit)
          break;  // while-loop breaks here if one v_kept is on boundary
        boundary_hit = true;
        nb = Neighbor(ff[0], edge.ord_in_faces[0], !ccw);
        continue;
      }

      nb.rotate(conn);
      v_kept_neighbor_edges.emplace_back(conn.face2edge[nb.f()][nb.j()]);
      v_kept_twins.emplace_back(conn.indices[nb.f()][nb.i()]);
      if (nb.f() == ff[0]) {
        // while-loop breaks here if there is no business of boundary
        assert(!boundary_hit);
        break;
      }
      // this face is deleted thus unnecessary to check fold-over if f == ff[0]
      v_kept_neighbors.emplace_back(nb);
    }

    // special case: there are 2 faces, 3 vertices in current component
    // this happens if input contains such component because
    // this edge collapse implementation avoid generating them
    if (v_del_twins.empty()) {
      assert(conn.indices[edge.faces[0]][edge.ord_in_faces[0]] ==
             conn.indices[edge.faces[1]][edge.ord_in_faces[1]]);
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

    assert(edge.ord_in_faces[0] != Edge::INVALID);
    const idx f = edge.faces[0];

    v_del_neighbors.reserve(ESTIMATE_VALENCE / 2);
    v_kept_neighbors.reserve(ESTIMATE_VALENCE / 2);
    v_del_twins.reserve(ESTIMATE_VALENCE / 2);
    v_kept_twins.reserve(ESTIMATE_VALENCE / 2);

    Neighbor nb(f, edge.ord_in_faces[0], ccw);
    while (true) {
      const auto& next_edge = conn.edge_of_face(nb.f(), nb.i());

      if (next_edge.on_boundary()) break;

      nb.rotate(conn);
      v_del_neighbors.emplace_back(nb);
      v_del_twins.emplace_back(conn.indices[nb.f()][nb.j()]);
    }

    nb = Neighbor(f, edge.ord_in_faces[0], !ccw);
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
    if (face_fold_over(vertices, conn.indices, nb_del, v_del, edge.center,
                       options.fold_over_angle_threshold))
      return false;
  for (const auto& nb_kept : v_kept_neighbors)
    if (face_fold_over(vertices, conn.indices, nb_kept, v_kept, edge.center,
                       options.fold_over_angle_threshold))
      return false;

  // check face quality
  for (const auto& nb_del : v_del_neighbors)
    if (extremely_elongated(vertices, conn.indices, nb_del, edge.center,
                            options.aspect_ratio_at_least))
      return false;
  for (const auto& nb_kept : v_del_neighbors)
    if (extremely_elongated(vertices, conn.indices, nb_kept, edge.center,
                            options.aspect_ratio_at_least))
      return false;

  // good to go
  return true;
}

bool edge_collapse(V& vertices, Internal::Connectivity& conn, Q& quadrics,
                   QEMHeap& heap, const idx ecol_target,
                   const SimplifyOptions& options) {
  auto& edge = conn.edges[ecol_target];

  // if non-boundary edge has two endpoints on boundary, we avoid collapsing it
  // because it is possible to produce non-manifold vertex
  if (edge.boundary_v == Edge::BOTH && !edge.on_boundary()) {
    heap.penalize(ecol_target);
    return false;
  }

  // collect neighboring faces into two containers
  std::vector<Neighbor> v_del_neighbors;
  std::vector<idx> v_kept_neighbor_edges;
  if (!scan_neighbors(vertices, conn, edge, v_del_neighbors,
                      v_kept_neighbor_edges, options)) {
    heap.penalize(ecol_target);
    return false;
  }

  // now that we are certain this edge is to be collapsed, remove it from heap
  heap.pop();

  const idx v_del = edge.vertices[edge.v_del_order()];
  const idx v_kept = edge.vertices[1 - edge.v_del_order()];

  if (!edge.on_boundary()) {
    // collapse a non-boundary edge and update its neighborhood

    const vec2i& ff = edge.faces;

    // update vertex position and quadric
    if (edge.boundary_v == Edge::NEITHER) vertices[v_kept] = edge.center;
    quadrics[v_kept] = edge.q;

    vec2i e_kept;  // global index of kept edge in two deleted faces
    for (auto i : {0, 1}) e_kept[i] = conn.edge_idx_across_from_v(ff[i], v_del);

    auto it = v_del_neighbors.begin();
    Edge* dirty_edge_ptr = &conn.edges[e_kept[0]];

    // first face to process: deleted face 0
    conn.indices[it->f()][it->center()] = v_kept;
    heap.erase(conn.face2edge[it->f()][it->j()]);
    conn.face2edge[it->f()][it->j()] = e_kept[0];
    const order f0_in_edge = dirty_edge_ptr->f_order(ff[0]);
    dirty_edge_ptr->faces[f0_in_edge] = it->f();
    dirty_edge_ptr->ord_in_faces[f0_in_edge] = it->j();

    // every face centered around deleted vertex
    for (++it; it != v_del_neighbors.end(); ++it) {
      conn.indices[it->f()][it->center()] = v_kept;
      dirty_edge_ptr = &conn.edge_of_face(it->f(), it->j());
      const order v_del_in_edge = dirty_edge_ptr->v_order(v_del);
      dirty_edge_ptr->vertices[v_del_in_edge] = v_kept;
      if (edge.boundary_v != Edge::NEITHER) {
        // not collapsing an edge entirely in the interior
        assert(edge.boundary_v == Edge::V0 || edge.boundary_v == Edge::V1);
        assert(dirty_edge_ptr->boundary_v != Edge::BOTH);
        if (dirty_edge_ptr->boundary_v == Edge::NEITHER)
          dirty_edge_ptr->boundary_v =
              static_cast<Edge::BOUNDARY_V>(v_del_in_edge);
        else
          dirty_edge_ptr->boundary_v = Edge::BOTH;
      }
      update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr,
                              options.fix_boundary);
    }

    // the deleted face 1
    --it;
    heap.erase(conn.face2edge[it->f()][it->i()]);
    conn.face2edge[it->f()][it->i()] = e_kept[1];
    dirty_edge_ptr = &conn.edges[e_kept[1]];
    const order f1_in_edge = dirty_edge_ptr->f_order(ff[1]);
    dirty_edge_ptr->faces[f1_in_edge] = it->f();
    dirty_edge_ptr->ord_in_faces[f1_in_edge] = it->i();

    // every edge centered around the kept vertex
    for (const idx e_dirty : v_kept_neighbor_edges)
      if (!(options.fix_boundary &&
            conn.edges[e_dirty].boundary_v == Edge::BOTH))
        update_error_and_center(vertices, quadrics, heap, &conn.edges[e_dirty],
                                options.fix_boundary);

  } else {
    // collapse a boundary edge and update its neighborhood

    const idx f = edge.faces[0];

    // update vertex position and quadric
    vertices[v_kept] = edge.center;
    quadrics[v_kept] = edge.q;

    const idx e_kept_idx = conn.edge_idx_across_from_v(f, v_del);
    Edge* dirty_edge_ptr = &conn.edges[e_kept_idx];

    // it is possible that v_del_neighbors is empty
    heap.erase(conn.edge_idx_across_from_v(f, v_kept));
    const order ford = dirty_edge_ptr->f_order(f);
    if (v_del_neighbors.empty()) {
      dirty_edge_ptr->ord_in_faces[ford] = Edge::INVALID;
      if (ford == 0) dirty_edge_ptr->swap_faces();
      assert(dirty_edge_ptr->ord_in_faces[0] != Edge::INVALID);
    } else {
      auto it = v_del_neighbors.begin();
      conn.face2edge[it->f()][it->j()] = e_kept_idx;
      dirty_edge_ptr->faces[ford] = it->f();
      dirty_edge_ptr->ord_in_faces[ford] = it->j();
    }

    // every face centered around deleted vertex
    for (const auto& nb : v_del_neighbors) {
      conn.indices[nb.f()][nb.center()] = v_kept;
      dirty_edge_ptr = &conn.edge_of_face(nb.f(), nb.i());
      const order v_del_in_edge = dirty_edge_ptr->v_order(v_del);
      dirty_edge_ptr->vertices[v_del_in_edge] = v_kept;
      assert(dirty_edge_ptr->boundary_v != Edge::NEITHER);
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
