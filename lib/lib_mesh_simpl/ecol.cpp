//
// Created by nickl on 1/8/19.
//

#include "ecol.h"
#include <algorithm>

namespace MeshSimpl {

namespace Internal {

static const size_t ESTIMATE_VALENCE = 8; // should cover more than 99.9% vertices

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
        edge.center = midpoint(vertices[edge.vertices[0]], vertices[edge.vertices[1]]);
        edge.error = q_error(edge.q, edge.center);
        for (const idx v : edge.vertices) {
            const double err = q_error(edge.q, vertices[v]);
            if (err < edge.error) {
                copy_vec3(vertices[v], edge.center);
                edge.error = err;
            }
        }
    }
}

void set_edge_error(const V& vertices, const Q& quadrics, Edge& edge) {
    assert(edge.boundary_v != Edge::BOTH);
    // error = v(Q1+Q2)v = vQv, v is new vertex position after edge-collapse
    const auto& vv = edge.vertices;
    edge.q = quadrics[vv[0]] + quadrics[vv[1]];

    if (edge.boundary_v == Edge::NONE)
        optimal_ecol_vertex_placement(vertices, edge);
    else {
        copy_vec3(vertices[vv[edge.boundary_v]], edge.center);
        edge.error = q_error(edge.q, edge.center);
    }
}

void compute_errors(const V& vertices, const Q& quadrics, E& edges) {
    // boundary edges will not be touched during simplification
    for (auto& edge : edges)
        if (edge.boundary_v != Edge::BOTH)
            set_edge_error(vertices, quadrics, edge);
}

bool face_fold_over(const V& vertices, const F& indices, const Neighbor& nb,
                    const idx v_move, const vec3d& move_to) {
    assert(v_move == indices[nb.f()][nb.center()]);
    const idx vi = indices[nb.f()][nb.i()];
    const idx vj = indices[nb.f()][nb.j()];
    const vec3d e0 = vertices[vj] - vertices[vi];
    const vec3d e1_prev = vertices[v_move] - vertices[vj];
    const vec3d e1_new = move_to - vertices[vj];
    vec3d normal_prev = cross(e0, e1_prev);
    double normal_prev_mag = magnitude(normal_prev);
    assert(normal_prev_mag != 0);
    vec3d normal_new = cross(e0, e1_new);
    double normal_new_mag = magnitude(normal_new);
    if (normal_new_mag == 0)
        return true;
    normal_prev /= normal_prev_mag;
    normal_new /= normal_new_mag;
    double cos = dot(normal_prev, normal_new);
    return cos < FOLD_OVER_COS_ANGLE;
}

void update_error_and_center(const V& vertices, const Q& quadrics, QEMHeap& heap,
                             Edge* const edge_ptr) {
    if (edge_ptr->boundary_v == Edge::BOTH) {
        heap.erase_by_ptr(edge_ptr);
    } else {
        const double error_prev = edge_ptr->error;
        set_edge_error(vertices, quadrics, *edge_ptr);
        heap.fix(edge_ptr, error_prev);
    }
}

bool scan_neighbors(const V& vertices, const F& indices, const E& edges,
                    const F2E& face2edge, const Edge& edge,
                    std::vector<Neighbor>& v_del_neighbors,
                    std::vector<idx>& v_kept_neighbor_edges) {
    const idx v_del = edge.vertices[choose_v_del(edge)];
    const idx v_kept = edge.vertices[1 - choose_v_del(edge)];
    std::vector<Neighbor> v_kept_neighbors;
    const vec2i& ff = edge.faces;
    std::vector<idx> v_del_twins, v_kept_twins;

    v_del_neighbors.reserve(ESTIMATE_VALENCE);
    v_kept_neighbors.reserve(ESTIMATE_VALENCE);
    v_del_twins.reserve(ESTIMATE_VALENCE);
    v_kept_twins.reserve(ESTIMATE_VALENCE);

    Neighbor nb{ff[0], edge.idx_in_face[0]};
    while (true) {
        nb.rotate(edges, face2edge);
        if (nb.f() == ff[1])
            break;
        v_del_neighbors.emplace_back(nb);
        v_del_twins.emplace_back(indices[nb.f()][nb.i()]);
    }

    bool boundary_hit = false; // useful in case v_kept is on boundary
    nb = {ff[1], edge.idx_in_face[1]};
    while (true) {
        const auto& next_edge = edges[face2edge[nb.f()][nb.i()]];

        if (next_edge.boundary_v == Edge::BOTH) {
            v_kept_twins.emplace_back(indices[nb.f()][nb.j()]);
            if (boundary_hit)
                break; // while-loop breaks here if one v_kept is on boundary
            boundary_hit = true;
            nb.counter_clockwise();
            nb = {ff[0], edge.idx_in_face[0]};
            continue;
        }

        nb.rotate(edges, face2edge);
        v_kept_neighbor_edges.emplace_back(face2edge[nb.f()][nb.j()]);
        v_kept_twins.emplace_back(indices[nb.f()][nb.i()]);
        if (nb.f() == ff[0])
            break; // while-loop breaks here if there is no business of boundary
        // this face is deleted thus unnecessary to check fold-over if f == ff[0]
        v_kept_neighbors.emplace_back(nb);
    }

    // check connectivity
    // Intersection check reference: en.cppreference.com/w/cpp/algorithm/set_intersection.html
    assert(!v_del_twins.empty());
    std::sort(v_del_twins.begin() + 1, v_del_twins.end());
    std::sort(v_kept_twins.begin(), v_kept_twins.end());
    for (auto i0 = v_del_twins.begin() + 1, i1 = v_kept_twins.begin();
         i0 != v_del_twins.end() && i1 != v_kept_twins.end();) {
        if (*i0 < *i1)
            ++i0;
        else if (*i0 > *i1)
            ++i1;
        else
            return false;
    }

    // check geometry
    for (const auto& nb_del : v_del_neighbors) {
        if (face_fold_over(vertices, indices, nb_del, v_del, edge.center))
            return false;
    }
    for (const auto& nb_kept : v_kept_neighbors) {
        if (face_fold_over(vertices, indices, nb_kept, v_kept, edge.center))
            return false;
    }

    // good to go
    return true;
}

bool edge_collapse(V& vertices, F& indices, E& edges, F2E& face2edge, Q& quadrics,
                   QEMHeap& heap, const idx ecol_target) {
    auto& edge = edges[ecol_target];
    const vec2i& ff = edge.faces;
    const idx v_del = edge.vertices[choose_v_del(edge)];
    const idx v_kept = edge.vertices[1 - choose_v_del(edge)];

    // to guarantee neighbors are traversed in clock-wise orientation,
    // swap faces (and idx_in_face) if necessary
    idx v_del_in_f0 = vi_in_face(indices, ff[0], v_del);
    if (edge.idx_in_face[0] != (v_del_in_f0 + 1) % 3) {
        std::swap(edge.idx_in_face[0], edge.idx_in_face[1]);
        std::swap(edge.faces[0], edge.faces[1]);
    }

    // collect neighboring faces into two containers
    std::vector<Neighbor> v_del_neighbors;
    std::vector<idx> v_kept_neighbor_edges;
    if (!scan_neighbors(vertices, indices, edges, face2edge, edge, v_del_neighbors,
                        v_kept_neighbor_edges)) {
        heap.penalize(ecol_target);
        return false;
    }

    // now that we are certain this edge is to be collapsed, remove it from heap
    heap.pop();
    // update vertex position and quadric
    if (edge.boundary_v == Edge::NONE)
        copy_vec3(edge.center, vertices[v_kept]);
    quadrics[v_kept] = edge.q;

    vec2i e_kept; // global index of kept edge in two deleted faces
    for (auto i : {0, 1})
        e_kept[i] = face2edge[ff[i]][vi_in_face(indices, ff[i], v_del)];

    auto it = v_del_neighbors.begin();
    Edge* dirty_edge_ptr = &edges[e_kept[0]];

    // first face to process: deleted face 0
    indices[it->f()][it->center()] = v_kept;
    heap.erase(face2edge[it->f()][it->j()]);
    face2edge[it->f()][it->j()] = e_kept[0];
    const idx ff0_in_edge = fi_in_edge(*dirty_edge_ptr, ff[0]);
    dirty_edge_ptr->faces[ff0_in_edge] = it->f();
    dirty_edge_ptr->idx_in_face[ff0_in_edge] = it->j();

    // every face centered around deleted vertex
    for (++it; it != v_del_neighbors.end(); ++it) {
        indices[it->f()][it->center()] = v_kept;
        dirty_edge_ptr = &edges[face2edge[it->f()][it->j()]];
        const idx v_del_in_edge = vi_in_edge(*dirty_edge_ptr, v_del);
        dirty_edge_ptr->vertices[v_del_in_edge] = v_kept;
        if (edge.boundary_v != Edge::NONE) {
            if (dirty_edge_ptr->boundary_v == Edge::NONE)
                dirty_edge_ptr->boundary_v = static_cast<Edge::BOUNDARY_V>(v_del_in_edge);
            else
                dirty_edge_ptr->boundary_v = Edge::BOTH;
        }
        update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr);
    }

    // the deleted face 1
    --it;
    heap.erase(face2edge[it->f()][it->i()]);
    face2edge[it->f()][it->i()] = e_kept[1];
    dirty_edge_ptr = &edges[e_kept[1]];
    const idx ff1_in_edge = fi_in_edge(*dirty_edge_ptr, ff[1]);
    dirty_edge_ptr->faces[ff1_in_edge] = it->f();
    dirty_edge_ptr->idx_in_face[ff1_in_edge] = it->i();

    // every edge centered around the kept vertex
    for (const idx e_dirty : v_kept_neighbor_edges)
        update_error_and_center(vertices, quadrics, heap, &edges[e_dirty]);

    return true;
}

} // namespace Internal

} // namespace MeshSimpl
