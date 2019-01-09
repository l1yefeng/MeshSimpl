//
// Created by nickl on 1/8/19.
//

#include "ecol.h"

namespace MeshSimpl
{

namespace Internal
{

void optimal_ecol_vertex_placement(const V& vertices, Edge& edge)
{
    const Quadric& q = edge.q;
    const vec3d b{q[6], q[7], q[8]};
    const double c = q[9];

    // computes the inverse of matrix A in quadric
    const double a_det = q[0]*(q[3]*q[5]-q[4]*q[4])
        -q[1]*(q[1]*q[5]-q[4]*q[2])
        +q[2]*(q[1]*q[4]-q[3]*q[2]);

    if (a_det != 0) {
        // invertible, find position yielding minimal error
        const double a_det_inv = 1.0/a_det;
        const std::array<double, 6> a_inv{(q[3]*q[5]-q[4]*q[4])*a_det_inv,
                                          (q[2]*q[4]-q[1]*q[5])*a_det_inv,
                                          (q[1]*q[4]-q[2]*q[3])*a_det_inv,
                                          (q[0]*q[5]-q[2]*q[2])*a_det_inv,
                                          (q[1]*q[2]-q[0]*q[4])*a_det_inv,
                                          (q[0]*q[3]-q[1]*q[1])*a_det_inv};
        edge.center = {-dot({a_inv[0], a_inv[1], a_inv[2]}, b),
                       -dot({a_inv[1], a_inv[3], a_inv[4]}, b),
                       -dot({a_inv[2], a_inv[4], a_inv[5]}, b)};
        edge.error = dot(b, edge.center)+c;
    }
    else {
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

void set_edge_error(const V& vertices, const std::vector<Quadric>& quadrics, Edge& edge)
{
    // error = v(Q1+Q2)v = vQv, v is new vertex position after edge-collapse
    const idx v0 = edge.vertices[0], v1 = edge.vertices[1];
    edge.q = quadrics[v0]+quadrics[v1];

    if (edge.boundary_v == BOUNDARY_V::NONE)
        optimal_ecol_vertex_placement(vertices, edge);
    else {
        copy_vec3(vertices[edge.boundary_v == BOUNDARY_V::V0 ? v0 : v1], edge.center);
        edge.error = q_error(edge.q, edge.center);
    }
}

void
init_edge_errors(const V& vertices, const std::vector<Quadric>& quadrics, std::vector<Edge>& edges)
{
    // boundary edges will not be touched during simplification
    for (auto& edge : edges)
        if (edge.boundary_v != BOUNDARY_V::BOTH)
            set_edge_error(vertices, quadrics, edge);

}

bool face_fold_over(const V& vertices, const idx v0, const idx v1, const idx v2_prev,
                    const vec3d& v2_new_pos)
{
    const vec3d e0 = vertices[v1]-vertices[v0];
    const vec3d e1_prev = vertices[v2_prev]-vertices[v1];
    const vec3d e1_new = v2_new_pos-vertices[v1];
    vec3d normal_prev = cross(e0, e1_prev);
    vec3d normal_new = cross(e0, e1_new);
    double normal_prev_mag = magnitude(normal_prev);
    if (normal_prev_mag == 0)
        return true;
    double normal_new_mag = magnitude(normal_new);
    if (normal_new_mag == 0)
        return true;
    double cos = dot(normal_prev /= normal_prev_mag, normal_new /= normal_new_mag);
    return cos < FOLD_OVER_COS_ANGLE;
}

void update_error_and_center(const V& vertices, const std::vector<Quadric>& quadrics, QEMHeap& heap,
                             Edge* const edge_ptr)
{
    edge_ptr->q = quadrics[edge_ptr->vertices[0]]+quadrics[edge_ptr->vertices[1]];
    const double error_prev = edge_ptr->error;
    optimal_ecol_vertex_placement(vertices, *edge_ptr);
    heap.fix(edge_ptr, edge_ptr->error > error_prev);
}

bool topology_preservation_check(const std::vector<vec3i>& face2edge,
                                 const std::queue<vec3i>& fve_queue_v_del,
                                 const std::queue<vec3i>& fve_queue_v_kept)
{
    vec3i fve;
    const idx& f = fve[0];
    const idx& v = fve[1];
    const idx& e = fve[2];
    fve = fve_queue_v_del.back();
    idx e_of_interest = face2edge[f][3-e-v];
    fve = fve_queue_v_kept.front();
    if (e_of_interest == face2edge[f][3-e-v])
        return false;
    fve = fve_queue_v_kept.back();
    e_of_interest = face2edge[f][3-e-v];
    fve = fve_queue_v_del.front();
    if (e_of_interest == face2edge[f][3-e-v])
        return false;
    return true;
}

bool collapse_interior_edge(V& vertices, F& indices, std::vector<Edge>& edges,
                            std::vector<vec3i>& face2edge, std::vector<Quadric>& quadrics,
                            std::vector<bool>& deleted_vertex, std::vector<bool>& deleted_face,
                            const std::function<void(vec3i&)>& iter_next, QEMHeap& heap,
                            const idx ecol_target)
{
    const auto& edge = edges[ecol_target];
    const vec2i& ff = edge.faces;
    // the choice of deletion does not matter
    idx v_kept = edge.vertices[0];
    idx v_del = edge.vertices[1];

    // local indexes of v_kept/v_del in two deleted faces
    vec2i v_kept_in_ff{vi_in_face(indices, ff[0], v_kept), vi_in_face(indices, ff[1], v_kept)};
    vec2i v_del_in_ff{vi_in_face(indices, ff[0], v_del), vi_in_face(indices, ff[1], v_del)};

    // global index of kept edge in two deleted faces
    vec2i e_kept{face2edge[ff[0]][v_del_in_ff[0]], face2edge[ff[1]][v_del_in_ff[1]]};
    vec3i fve;
    const idx& f = fve[0];
    const idx& v = fve[1];
    const idx& e = fve[2]; // they always refer to fve; updated on each `fve = ...'

    // starting point of fve iterations at two deleted faces
    const vec3i fve_begin_v_del{ff[0], v_kept_in_ff[0], edge.idx_in_face[0]};
    const vec3i fve_begin_v_kept{ff[1], v_del_in_ff[1], edge.idx_in_face[1]};

    // test run iterating faces: need to increase error and abort if fold-over is identified
    std::queue<vec3i> fve_queue_v_del, fve_queue_v_kept;
    for (iter_next(fve = fve_begin_v_del); f != ff[1]; iter_next(fve)) {
        fve_queue_v_del.push(fve);
        if (Internal::face_fold_over(vertices, indices[f][e], indices[f][v], v_del, edge.center))
            return false;
    }
    for (iter_next(fve = fve_begin_v_kept); f != ff[0]; iter_next(fve)) {
        fve_queue_v_kept.push(fve);
        if (Internal::face_fold_over(vertices, indices[f][e], indices[f][v], v_kept, edge.center))
            return false;
    }

    // check if this operation will damage topology (create 3-manifold)
    if (!topology_preservation_check(face2edge, fve_queue_v_del, fve_queue_v_kept))
        return false;

    // mark 2 faces, 1 vertex, and 1 edge as deleted;
    // 2 more edges will be marked in the loop later
    deleted_face[ff[0]] = true;
    deleted_face[ff[1]] = true;
    deleted_vertex[v_del] = true;
    heap.pop();
    // update vertex position and quadric
    copy_vec3(edge.center, vertices[v_kept]);
    quadrics[v_kept] = edge.q;

    fve = fve_queue_v_del.front();
    Edge* tgt_edge = &edges[e_kept[0]];
    // first face to process: deleted face 0
    indices[f][3-v-e] = v_kept;
    heap.erase(face2edge[f][e]);
    face2edge[f][e] = e_kept[0];
    const idx ff0_in_edge = fi_in_edge(*tgt_edge, ff[0]);
    assert(ff[0] != f);
    tgt_edge->faces[ff0_in_edge] = f;
    tgt_edge->idx_in_face[ff0_in_edge] = e;
    // every face centered around deleted vertex
    fve_queue_v_del.pop();
    for (; !fve_queue_v_del.empty(); fve_queue_v_del.pop()) {
        fve = fve_queue_v_del.front();
        indices[f][3-v-e] = v_kept;
        tgt_edge = &edges[face2edge[f][e]];
        tgt_edge->vertices = {indices[f][v], v_kept};
        update_error_and_center(vertices, quadrics, heap, tgt_edge);
    }
    // the deleted face 1
    heap.erase(face2edge[f][v]);
    face2edge[f][v] = e_kept[1];
    tgt_edge = &edges[e_kept[1]];
    const idx ff1_in_edge = fi_in_edge(*tgt_edge, ff[1]);
    assert(ff[1] != f);
    tgt_edge->faces[ff1_in_edge] = f;
    tgt_edge->idx_in_face[ff1_in_edge] = v;

    // every face centered around the kept vertex
    for (; !fve_queue_v_kept.empty(); fve_queue_v_kept.pop()) {
        fve = fve_queue_v_kept.front();
        tgt_edge = &edges[face2edge[f][e]];
        update_error_and_center(vertices, quadrics, heap, tgt_edge);
    }
    assert(face2edge[f][v] == e_kept[0]);
    tgt_edge = &edges[face2edge[f][v]];
    update_error_and_center(vertices, quadrics, heap, tgt_edge);

    return true;
}

}

}