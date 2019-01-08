#ifndef MESH_SIMPL_LIBRARY_H
#define MESH_SIMPL_LIBRARY_H

#include "types.h"
#include "util.h"
#include "qem_heap.h"
#include <vector>
#include <array>
#include <set>
#include <queue>
#include <limits>
#include <functional>

namespace MeshSimpl
{

namespace Internal
{

static const size_t MIN_NR_VERTICES = 4;

// Compute quadrics Q for every vertex. If not weighted by area then it's uniform
std::vector<Quadric>
compute_quadrics(const V& vertices, const F& indices, bool weight_by_area)
{
    // quadrics are initialized with all zeros
    std::vector<Quadric> quadrics(vertices.size());

    for (const auto& face : indices) {
        // calculate the plane of this face (n and d: n'v+d=0 defines the plane)
        const vec3d edge1 = vertices[face[1]]-vertices[face[0]];
        const vec3d edge2 = vertices[face[2]]-vertices[face[0]];
        vec3d normal = cross(edge1, edge2);
        // |normal| = area, used for normalization and weighting quadrics
        double area = magnitude(normal);
        normal /= area;
        // d = -n*v0
        double d = -dot(normal, vertices[face[0]]);

        // calculate quadric Q = (A, b, c) = (nn', dn, d*d)
        Quadric q = make_quadric(normal, d);

        if (weight_by_area)
            q *= area;

        for (const auto& v : face)
            quadrics[v] += q;
    }

    return quadrics;
}

// Returns face2edge and edges. face2edge is NFx3 with each value indexing a unique edge in edges.
// This method does not initialize members optimal_pos and error in struct Edge.
std::pair<std::vector<Edge>, std::vector<vec3i>>
edge_topology(const F& indices, const size_t vertex_cnt)
{
    const auto edge_cmp = [](const Edge& a, const Edge& b)->bool {
        if (a.vertices[0] < b.vertices[0])
            return true;
        if (a.vertices[0] > b.vertices[0])
            return false;
        return a.vertices[1] < b.vertices[1];
    };
    std::set<Edge, decltype(edge_cmp)> edge_set(edge_cmp);

    // insert all edges into edge_set and find out if it is on boundary
    for (idx f = 0; f < indices.size(); ++f) {
        const auto& face = indices[f];
        for (idx i = 0; i < 3; ++i) {
            const idx j = (i+1)%3; // face[i] and face[j] forms this edge
            const idx k = (i+2)%3; // k is index of the other vertex and of this edge local to face
            idx v0 = face[i], v1 = face[j];
            if (v0 > v1)
                std::swap(v0, v1);
            auto it_and_inserted = edge_set.insert({{v0, v1}, {f}, {k}, BOUNDARY_V::BOTH});
            auto it = it_and_inserted.first;
            if (!it_and_inserted.second) {
                // edge was not inserted because it is already there
                // modifying through immutable iterator only if do not affect order
                *const_cast<BOUNDARY_V*>(&it->boundary_v) = BOUNDARY_V::NONE;
                *const_cast<idx*>(&it->faces[1]) = f;
                *const_cast<idx*>(&it->idx_in_face[1]) = k;
            }
        }
    }

    // convert edges from set to vector
    std::vector<Edge> edges(edge_set.begin(), edge_set.end());
    std::vector<bool> vertex_on_boundary(vertex_cnt, false);
    std::vector<vec3i> face2edge(indices.size());

    for (idx i = 0; i < edges.size(); ++i) {
        const auto& edge = edges[i];
        // identify boundary vertices
        if (edge.boundary_v == BOUNDARY_V::BOTH) {
            vertex_on_boundary[edge.vertices[0]] = true;
            vertex_on_boundary[edge.vertices[1]] = true;
        }

        // populate face2edge references
        for (auto j : {0, 1})
            face2edge[edge.faces[j]][edge.idx_in_face[j]] = i;
    }

    // non-boundary edges may have one vertex on boundary, find them in this loop
    for (auto& edge : edges) {
        if (edge.boundary_v == BOUNDARY_V::BOTH)
            continue;
        if (vertex_on_boundary[edge.vertices[0]] && vertex_on_boundary[edge.vertices[1]])
            edge.boundary_v = BOUNDARY_V::BOTH;
        else if (vertex_on_boundary[edge.vertices[0]])
            edge.boundary_v = BOUNDARY_V::V0;
        else if (vertex_on_boundary[edge.vertices[1]])
            edge.boundary_v = BOUNDARY_V::V1;
        // else totally within the boundary
    }

    return {edges, face2edge};
}

// Set error and center of edge by choosing a position to collapse into
void ecol_vertex_placement(const V& vertices, Edge& edge)
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
                copy_vertex_position(vertices[v], edge.center);
                edge.error = err;
            }
        }
    }
}

// Compute quadric error for every edge at initialization.
// Returns the range of error, can be used to decide fold-over penalty factor
double init_edge_errors(const V& vertices,
                        const std::vector<Quadric>& quadrics,
                        std::vector<Edge>& edges)
{
    double min_error = std::numeric_limits<double>::max();
    double max_error = std::numeric_limits<double>::lowest();
    for (auto& edge : edges) {
        // boundary edges will not be touched during simplification
        if (edge.boundary_v == BOUNDARY_V::BOTH)
            continue;

        // error = v(Q1+Q2)v = vQv, v is new vertex position after edge-collapse
        const idx v0 = edge.vertices[0], v1 = edge.vertices[1];
        edge.q = quadrics[v0]+quadrics[v1];

        if (edge.boundary_v == BOUNDARY_V::NONE)
            ecol_vertex_placement(vertices, edge);
        else {
            copy_vertex_position(vertices[edge.boundary_v == BOUNDARY_V::V0 ? v0 : v1],
                                 edge.center);
            edge.error = q_error(edge.q, edge.center);
        }
        min_error = std::min(min_error, edge.error);
        max_error = std::max(max_error, edge.error);
    }

    return max_error-min_error;
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
    return cos < -0.97562931279;
}

// Replace a vertex in edge, and update other members then fix priority in heap
void update_error_and_center(const V& vertices, const std::vector<Quadric>& quadrics,
                             QEMHeap& heap, Edge* const edge_ptr)
{
    edge_ptr->q = quadrics[edge_ptr->vertices[0]]+quadrics[edge_ptr->vertices[1]];
    const double error_prev = edge_ptr->error;
    ecol_vertex_placement(vertices, *edge_ptr);
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

// Remove the vertices and indices that are marked deleted, and reduce the vector size
void compact_data(const std::vector<bool>& deleted_vertex, const std::vector<bool>& deleted_face,
                  V& vertices, F& indices)
{
    std::vector<std::array<std::vector<idx>, 3>> vertex2face(vertices.size());

    // get rid of all deleted faces
    for (size_t lo = 0, hi = indices.size()-1;; ++lo, --hi) {
        while (!deleted_face[lo] && lo <= hi)
            ++lo;
        while (deleted_face[hi] && lo < hi)
            --hi;
        if (lo >= hi) {
            indices.resize(lo);
            break;
        }
        std::swap(indices[lo], indices[hi]);
    }

    // create mapping from v to f
    for (idx f = 0; f < indices.size(); ++f)
        for (idx i = 0; i < 3; ++i)
            vertex2face[indices[f][i]][i].push_back(static_cast<idx>(f));

    // get rid of deleted vertices and keep the mapping valid
    for (size_t lo = 0, hi = vertices.size()-1;; ++lo, --hi) {
        while (!deleted_vertex[lo] && lo <= hi)
            ++lo;
        while (deleted_vertex[hi] && lo < hi)
            --hi;
        if (lo >= hi) {
            vertices.resize(lo);
            vertex2face.resize(lo);
            break;
        }
        std::swap(vertex2face[lo], vertex2face[hi]);
        std::swap(vertices[lo], vertices[hi]);
    }

    for (idx v = 0; v < vertex2face.size(); ++v)
        for (idx i = 0; i < 3; ++i)
            for (const idx f : vertex2face[v][i])
                indices[f][i] = static_cast<idx>(v);
}

inline idx vi_in_face(const F& indices, const idx f, const idx v)
{
    if (indices[f][0] == v)
        return 0;
    if (indices[f][1] == v)
        return 1;
    assert(indices[f][2] == v);
    return 2;
}

inline idx fi_in_edge(const Edge& edge, const idx f)
{
    if (edge.faces[0] == f)
        return 0;
    assert(edge.faces[1] == f);
    return 1;
};

inline idx vi_in_edge(const Edge& edge, const idx v)
{
    if (edge.vertices[0] == v)
        return 0;
    assert(edge.vertices[1] == v);
    return 1;
}

// Returns true if edge is collapsed
bool collapse_interior_edge(V& vertices, F& indices,
                            std::vector<Edge>& edges,
                            std::vector<vec3i>& face2edge,
                            std::vector<Quadric>& quadrics,
                            std::vector<bool>& deleted_vertex,
                            std::vector<bool>& deleted_face,
                            const std::function<void(vec3i&)>& iter_next,
                            QEMHeap& heap, const idx ecol_target)
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
    copy_vertex_position(edge.center, vertices[v_kept]);
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

} // namespace MeshSimpl::Internal

// Mesh simplification main method. Simplify given mesh until remaining number of vertices/faces
// is (1-strength) of the original. Returns output vertices and indices as in inputs.
// TODO: Ignoring the 4th and the following values (if exist) in vertices.
std::pair<V, F> simplify(const V& vertices, const F& indices, float strength)
{
    const auto nv_target = std::max(static_cast<size_t>(std::lround((1-strength)*vertices.size())),
                                    Internal::MIN_NR_VERTICES);
    auto out_vertices = vertices;
    auto out_indices = indices;

    if (nv_target == vertices.size())
        return {out_vertices, out_indices};

    auto quadrics = Internal::compute_quadrics(vertices, indices, true);

    auto edge_topo = Internal::edge_topology(indices, vertices.size());
    std::vector<Internal::Edge>& edges = edge_topo.first;
    std::vector<vec3i>& face2edge = edge_topo.second;

    const double error_penalty_factor = Internal::init_edge_errors(vertices, quadrics, edges)/100;

    // create priority queue on quadric error
    Internal::QEMHeap heap(edges);

    std::vector<bool> deleted_vertex(vertices.size(), false);
    std::vector<bool> deleted_face(indices.size(), false);

    // returns f,v,e of next neighbor face centered around v_center;
    // indexes v and e are local in face (0,1,2) in both input and output.
    const auto iter_next = [&](vec3i& fve)->void {
        const idx f = fve[0], v = fve[1], e = fve[2];
        assert(v != e);
        const auto& edge = edges[face2edge[f][v]];
        const idx f_idx_to_edge = fi_in_edge(edge, f);
        const idx of = edge.faces[1-f_idx_to_edge];
        assert(f != of);
        const idx ov_global = out_indices[f][e];
        const idx ov = Internal::vi_in_face(out_indices, of, ov_global);
        assert(face2edge[of][ov] != face2edge[f][e]);
        const idx oe = edge.idx_in_face[1-f_idx_to_edge];
        fve = {of, ov, oe};
    };

    for (auto nv = vertices.size(); !heap.empty() && nv > nv_target;) {
        // collapse an edge to remove 1 vertex and 2 faces in each iteration
        const idx e_collapsed = heap.top();

        const auto& edge = edges[e_collapsed];
        assert(edge.vertices[0] != edge.vertices[1]);
        assert(edge.faces[0] != edge.faces[1]);

        if (edge.boundary_v == Internal::BOUNDARY_V::NONE) {
            if (collapse_interior_edge(out_vertices, out_indices, edges, face2edge, quadrics,
                                       deleted_vertex, deleted_face,
                                       iter_next, heap, e_collapsed)) {
                --nv;
            }
            else {
                // cause of failure is fold-over, let's penalize it and give it a second chance
                heap.penalize(e_collapsed);
            }
        }
        else {
            // TODO: handle boundary cases
            assert(false);
        }
    }

    // edges are pointless from this point on, but need to fix vertices and indices
    Internal::compact_data(deleted_vertex, deleted_face, out_vertices, out_indices);
    return {out_vertices, out_indices};
}

} // namespace MeshSimpl

#endif // MESH_SIMPL_LIBRARY_H
