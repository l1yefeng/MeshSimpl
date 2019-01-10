//
// Created by nickl on 1/8/19.
//

#include "ecol.h"
#include <algorithm>

namespace MeshSimpl {

namespace Internal {

void optimal_ecol_vertex_placement(const V& vertices, Edge& edge) {
    const Quadric& q = edge.q;
    const vec3d b{q[6], q[7], q[8]};
    const double c = q[9];

    // computes the inverse of matrix A in quadric
    const double a_det = q[0] * (q[3] * q[5] - q[4] * q[4]) - q[1] * (q[1] * q[5] - q[4] * q[2]) +
                         q[2] * (q[1] * q[4] - q[3] * q[2]);

    if (a_det != 0) {
        // invertible, find position yielding minimal error
        const double a_det_inv = 1.0 / a_det;
        const std::array<double, 6> a_inv{
            (q[3] * q[5] - q[4] * q[4]) * a_det_inv, (q[2] * q[4] - q[1] * q[5]) * a_det_inv,
            (q[1] * q[4] - q[2] * q[3]) * a_det_inv, (q[0] * q[5] - q[2] * q[2]) * a_det_inv,
            (q[1] * q[2] - q[0] * q[4]) * a_det_inv, (q[0] * q[3] - q[1] * q[1]) * a_det_inv};
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
    // error = v(Q1+Q2)v = vQv, v is new vertex position after edge-collapse
    const idx v0 = edge.vertices[0], v1 = edge.vertices[1];
    edge.q = quadrics[v0] + quadrics[v1];

    if (edge.boundary_v == BOUNDARY_V::NONE)
        optimal_ecol_vertex_placement(vertices, edge);
    else {
        copy_vec3(vertices[edge.boundary_v == BOUNDARY_V::V0 ? v0 : v1], edge.center);
        edge.error = q_error(edge.q, edge.center);
    }
}

void compute_errors(const V& vertices, const Q& quadrics, E& edges) {
    // boundary edges will not be touched during simplification
    for (auto& edge : edges)
        if (edge.boundary_v != BOUNDARY_V::BOTH)
            set_edge_error(vertices, quadrics, edge);
}

void recompute_errors(const V& vertices, const Q& quadrics, E& edges,
                      std::vector<idx>::const_iterator edge_idx_begin,
                      std::vector<idx>::const_iterator edge_idx_end) {
    for (auto it = edge_idx_begin; it != edge_idx_end; ++it)
        set_edge_error(vertices, quadrics, edges[*it]);
}

bool face_fold_over(const V& vertices, const idx v0, const idx v1, const idx v2_prev,
                    const vec3d& v2_new_pos) {
    const vec3d e0 = vertices[v1] - vertices[v0];
    const vec3d e1_prev = vertices[v2_prev] - vertices[v1];
    const vec3d e1_new = v2_new_pos - vertices[v1];
    vec3d normal_prev = cross(e0, e1_prev);
    double normal_prev_mag = magnitude(normal_prev);
    assert(normal_prev_mag != 0);
    vec3d normal_new = cross(e0, e1_new);
    double normal_new_mag = magnitude(normal_new);
    if (normal_new_mag == 0)
        return true;
    double cos = dot(normal_prev /= normal_prev_mag, normal_new /= normal_new_mag);
    return cos < FOLD_OVER_COS_ANGLE;
}

void update_error_and_center(const V& vertices, const Q& quadrics, QEMHeap& heap,
                             Edge* const edge_ptr) {
    edge_ptr->q = quadrics[edge_ptr->vertices[0]] + quadrics[edge_ptr->vertices[1]];
    const double error_prev = edge_ptr->error;
    optimal_ecol_vertex_placement(vertices, *edge_ptr);
    heap.fix(edge_ptr, edge_ptr->error > error_prev);
}

bool scan_neighbors(const V& vertices, const F& indices, const E& edges, const F2E& face2edge,
                    const Edge& edge, const idx v_del, const idx v_kept,
                    std::vector<vec3i>& star_v_del, std::vector<vec3i>& star_v_kept) {
    const vec2i& ff = edge.faces;
    std::vector<idx> v0_twins, v1_twins;

    idx f = ff[0];
    idx v = vi_in_face(indices, ff[0], v_kept);
    idx e = edge.idx_in_face[0];
    std::vector<vec3i>* star_ptr = &star_v_del;
    while (true) {
        assert(v != e);
        const auto& curr_edge = edges[face2edge[f][v]];
        const idx f_idx_to_edge = fi_in_edge(curr_edge, f);
        const idx of = curr_edge.faces[1 - f_idx_to_edge];
        assert(f != of);
        const idx ov_global = indices[f][e];
        const idx ov = Internal::vi_in_face(indices, of, ov_global);
        assert(face2edge[of][ov] != face2edge[f][e]);
        const idx oe = curr_edge.idx_in_face[1 - f_idx_to_edge];
        f = of;
        v = ov;
        e = oe;

        if (f == ff[1]) {
            f = ff[1];
            v = vi_in_face(indices, ff[1], v_del);
            e = edge.idx_in_face[1];
            star_ptr = &star_v_kept;
            std::swap(v0_twins, v1_twins);
            continue;
        } else if (f == ff[0]) {
            break;
        }
        star_ptr->push_back({f, v, e});
        v0_twins.push_back(ov_global);
    }

    // check connectivity
    std::sort(v0_twins.begin(), v0_twins.end());
    std::sort(v1_twins.begin(), v1_twins.end());
    for (auto i0 = v0_twins.begin(), i1 = v1_twins.begin();
         i0 != v0_twins.end() && i1 != v1_twins.end();) {
        if (*i0 < *i1)
            ++i0;
        else if (*i0 > *i1)
            ++i1;
        else
            return false;
    }

    // check geometry
    for (const auto& fve : star_v_del) {
        if (face_fold_over(vertices, indices[f][fve[2]], indices[f][fve[1]], v_del, edge.center))
            return false;
    }
    for (const auto& fve : star_v_kept) {
        if (face_fold_over(vertices, indices[f][fve[2]], indices[f][fve[1]], v_kept, edge.center))
            return false;
    }

    // good to go
    return true;
}

bool collapse_interior_edge(V& vertices, F& indices, E& edges, F2E& face2edge, Q& quadrics,
                            std::vector<bool>& deleted_vertex, std::vector<bool>& deleted_face,
                            QEMHeap& heap, const idx ecol_target) {
    const auto& edge = edges[ecol_target];
    const vec2i& ff = edge.faces;
    // the choice of deletion does not matter
    const idx v_kept = edge.vertices[0];
    const idx v_del = edge.vertices[1];

    // collect neighboring faces into fifo queues and meanwhile detect possible fold-over problems
    std::vector<vec3i> star_v_del, star_v_kept;
    if (!scan_neighbors(vertices, indices, edges, face2edge, edge, v_del, v_kept, star_v_del,
                        star_v_kept)) {
        heap.penalize(ecol_target);
        return false;
    }

    // mark 2 faces and 1 vertex as deleted
    deleted_face[ff[0]] = true;
    deleted_face[ff[1]] = true;
    deleted_vertex[v_del] = true;
    // now that we are certain this edge is to be collapsed, remove it from heap
    heap.pop();
    // update vertex position and quadric
    copy_vec3(edge.center, vertices[v_kept]);
    quadrics[v_kept] = edge.q;

    vec2i e_kept; // global index of kept edge in two deleted faces
    for (auto i : {0, 1})
        e_kept[i] = face2edge[ff[i]][vi_in_face(indices, ff[i], v_del)];
    vec3i fve;
    const idx& f = fve[0];
    const idx& v = fve[1];
    const idx& e = fve[2]; // they always refer to fve; updated on each `fve = ...'

    fve = star_v_del[0];
    Edge* dirty_edge_ptr = &edges[e_kept[0]];
    // first face to process: deleted face 0
    indices[f][fve_center(fve)] = v_kept;
    heap.erase(face2edge[f][e]);
    face2edge[f][e] = e_kept[0];
    const idx ff0_in_edge = fi_in_edge(*dirty_edge_ptr, ff[0]);
    assert(ff[0] != f);
    dirty_edge_ptr->faces[ff0_in_edge] = f;
    dirty_edge_ptr->idx_in_face[ff0_in_edge] = e;
    // every face centered around deleted vertex
    for (auto it = star_v_del.begin() + 1; it != star_v_del.end(); ++it) {
        fve = *it;
        indices[f][fve_center(fve)] = v_kept;
        dirty_edge_ptr = &edges[face2edge[f][e]];
        dirty_edge_ptr->vertices = {indices[f][v], v_kept};
        update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr);
    }
    // the deleted face 1
    heap.erase(face2edge[f][v]);
    face2edge[f][v] = e_kept[1];
    dirty_edge_ptr = &edges[e_kept[1]];
    const idx ff1_in_edge = fi_in_edge(*dirty_edge_ptr, ff[1]);
    assert(ff[1] != f);
    dirty_edge_ptr->faces[ff1_in_edge] = f;
    dirty_edge_ptr->idx_in_face[ff1_in_edge] = v;

    // every face centered around the kept vertex
    for (auto& _fve : star_v_kept) {
        fve = _fve;
        dirty_edge_ptr = &edges[face2edge[f][e]];
        update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr);
    }
    assert(face2edge[f][v] == e_kept[0]);
    dirty_edge_ptr = &edges[face2edge[f][v]];
    update_error_and_center(vertices, quadrics, heap, dirty_edge_ptr);

    return true;
}

} // namespace Internal

} // namespace MeshSimpl
