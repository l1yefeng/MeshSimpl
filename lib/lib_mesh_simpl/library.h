#ifndef MESH_SIMPL_LIBRARY_H
#define MESH_SIMPL_LIBRARY_H

#include "util.h"
#include <vector>
#include <array>
#include <set>
#include <queue>
#include <algorithm>

namespace MeshSimpl
{

namespace Internal
{

static const double FOLD_OVER_PENALTY = 1.5;
static const size_t MIN_NR_VERTICES = 4;

enum BOUNDARY_V { NONE, BOTH, V0, V1 };

// Edge defines the struct of an edge
struct Edge
{
    vec2i vertices;        // index of two end vertices, unique, v0 < v1
    vec2i faces;           // index of two incident faces
    vec2i idx_in_face;     // index of this (0, 1, 2) in faces
    BOUNDARY_V boundary_v; // vertices on boundary
    vec3d center;          // where this edge collapse into
    double error;          // quadric error value, undefined if on boundary
    Quadric q;             // sum of quadrics of two vertices
    bool dirty;            // true if this edge has wrong position in heap
};

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
        for (idx j = 0; j < 2; ++j)
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

// Compute quadric error for every edge at initialization.
void init_edge_errors(const V& vertices,
                      const std::vector<Quadric>& quadrics,
                      std::vector<Edge>& edges)
{
    for (auto& edge : edges) {
        // boundary edges will not be touched during simplification
        if (edge.boundary_v == BOUNDARY_V::BOTH)
            continue;

        // error = v(Q1+Q2)v = vQv, v is new vertex position after edge-collapse
        const auto& vv = edge.vertices;
        Quadric q = quadrics[vv[0]]+quadrics[vv[1]];

        if (edge.boundary_v == BOUNDARY_V::V0)
            std::copy(vertices[vv[0]].begin(), vertices[vv[0]].end(), edge.center.begin());
        else if (edge.boundary_v == BOUNDARY_V::V1)
            std::copy(vertices[vv[1]].begin(), vertices[vv[1]].end(), edge.center.begin());
        else // edge.optimal_pos == V_ON_BOUNDARY::NONE
            edge.center = optimal_v_pos(q);

        edge.error = q_error(q, edge.center, edge.boundary_v == BOUNDARY_V::NONE);
        edge.q = q;
    }
}

// Returns a min-heap consists of all non-boundary edge indexes.
// Using auto keyword due to the lambda compare function -- requires C++14
auto build_min_heap(std::vector<Edge>& edges)
{
    const auto cmp = [&](idx a, idx b)->bool { return edges[a].error > edges[b].error; };
    std::priority_queue<idx, std::vector<idx>, decltype(cmp)> heap(cmp);
    for (idx i = 0; i < edges.size(); ++i)
        if (edges[i].boundary_v != BOUNDARY_V::BOTH)
            heap.push(i);
    return heap;
}

bool face_fold_over(const V& vertices, const idx v0, const idx v1, const idx v2_prev,
                    const vec3d& v2_new_pos)
{
    const vec3d e0 = vertices[v1]-vertices[v0];
    const vec3d e1_prev = vertices[v2_prev]-vertices[v1];
    const vec3d e1_new = v2_new_pos-vertices[v1];
    const vec3d normal_prev = cross(e0, e1_prev);
    const vec3d normal_new = cross(e0, e1_new);
    return dot(normal_prev, normal_new) < 0;
}

// Replace a vertex in edge, and update other members then mark dirty.
// FIXME: only for non-boundary cases for the time being
void update_edge(Edge& edge, const Quadric& q0, const Quadric& q1)
{
    edge.q = q0+q1;
    edge.center = optimal_v_pos(edge.q);
    edge.error = q_error(edge.q, edge.center, true);
    edge.dirty = true;
}

// Replace a vertex in edge, and update other members then mark dirty.
// FIXME: only for non-boundary cases for the time being
void update_edge(Edge& edge, const Quadric& q0, const Quadric& q1, const idx v0, const idx v1)
{
    edge.vertices[0] = v0;
    edge.vertices[1] = v1;
    update_edge(edge, q0, q1);
}

// Remove the vertices and indices that are marked deleted, and reduce the vector size
void compact_data(V& vertices, F& indices,
                  const std::vector<bool>& vertex_deleted,
                  const std::vector<bool>& face_deleted)
{
    std::vector<std::array<std::vector<idx>, 3>> vertex2face(vertices.size());

    // get rid of all deleted faces
    for (size_t lo = 0, hi = indices.size()-1;; ++lo, --hi) {
        while (!face_deleted[lo] && lo <= hi)
            ++lo;
        while (face_deleted[hi] && lo < hi)
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
        while (!vertex_deleted[lo] && lo <= hi)
            ++lo;
        while (vertex_deleted[hi] && lo < hi)
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

    auto quadrics = Internal::compute_quadrics(vertices, indices, true);

    auto edge_topo = Internal::edge_topology(indices, vertices.size());
    std::vector<Internal::Edge>& edges = edge_topo.first;
    std::vector<vec3i>& face2edge = edge_topo.second;

    Internal::init_edge_errors(vertices, quadrics, edges);

    // create priority queue on quadric error
    auto heap = Internal::build_min_heap(edges);
    std::vector<bool> vertex_deleted(vertices.size(), false);
    std::vector<bool> face_deleted(indices.size(), false);
    std::vector<bool> edge_deleted(edges.size(), false);

    // returns the position of vertex in face. a convenient method but make sure v EXISTS in face
    const auto vi_in_face = [&](const idx f, const idx v)->idx {
        const auto& row = out_indices[f];
        auto i = std::distance(row.begin(), std::find(row.begin(), row.end(), v));
        assert(i < 3);
        return static_cast<idx>(i);
    };

    const auto fi_in_edge = [](const Internal::Edge& edge, const idx f)->idx {
        return edge.faces[0] == f ? 0 : 1;
    };

    // returns f,v,e of next neighbor face centered around v_center;
    // indexes v and e are local in face (0,1,2) in both input and output.
    const auto iter_next = [&](const vec3i& fve)->vec3i {
        const idx f = fve[0], v = fve[1], e = fve[2];
        assert(v != e);
        const auto& edge = edges[face2edge[f][v]];
        const idx f_idx_to_edge = edge.faces[0] == f ? 0 : 1;
        const idx of = edge.faces[1-f_idx_to_edge];
        const idx oe = edge.idx_in_face[1-f_idx_to_edge];
        const idx ov_global = out_indices[f][e];
        const idx ov = vi_in_face(of, ov_global);
        return {of, ov, oe};
    };

    for (auto nv = vertices.size(); !heap.empty() && nv > nv_target;) {
        // collapse an edge to remove 1 vertex and 2 faces in each iteration
        const idx e_collapsed = heap.top();
        heap.pop();

        if (edge_deleted[e_collapsed])
            continue;

        const auto& edge = edges[e_collapsed];
        if (edge.dirty) {
            // this edge was modified and is having an error no less than before
            edges[e_collapsed].dirty = false;
            heap.push(e_collapsed);
            continue;
        }

        const auto& vv = edge.vertices; // one is kept and another will be deleted
        const auto& ff = edge.faces;    // two faces that will be deleted

        if (edge.boundary_v == Internal::BOUNDARY_V::NONE) {
            // the choice of deletion does not matter
            const idx v_kept = vv[0];
            const idx v_del = vv[1];

            // local indexes of v_kept/v_del in two deleted faces
            const vec2i v_kept_in_ff{vi_in_face(ff[0], v_kept), vi_in_face(ff[1], v_kept)};
            const vec2i v_del_in_ff{vi_in_face(ff[0], v_del), vi_in_face(ff[1], v_del)};

            // starting point of fve iterations at two deleted faces
            const vec3i fve_begin_v_del{ff[0], v_kept_in_ff[0], edge.idx_in_face[0]};
            const vec3i fve_begin_v_kept{ff[1], v_del_in_ff[1], edge.idx_in_face[1]};

            // test run iterating faces: need to increase error and abort if fold-over is identified
            bool danger_of_fold_over = false;
            std::queue<vec3i> fve_queue_v_del, fve_queue_v_kept;
            vec3i fve{};
            const idx& f = fve[0];
            const idx& v = fve[1];
            const idx& e = fve[2]; // they always refer to fve; updated on each `fve = ...'
            for (fve = iter_next(fve_begin_v_del); f != ff[1]; fve = iter_next(fve)) {
                fve_queue_v_del.push(fve);
                if (Internal::face_fold_over(out_vertices, out_indices[f][e], out_indices[f][v],
                                             v_del, edge.center)) {
                    danger_of_fold_over = true;
                    break;
                }
            }
            if (!danger_of_fold_over) {
                for (fve = iter_next(fve_begin_v_kept); f != ff[0]; fve = iter_next(fve)) {
                    fve_queue_v_kept.push(fve);
                    if (Internal::face_fold_over(out_vertices, out_indices[f][e], out_indices[f][v],
                                                 v_kept, edge.center)) {
                        danger_of_fold_over = true;
                        break;
                    }
                }
            }
            if (danger_of_fold_over) {
                // penalize this edge by increasing its error and then push back to min-heap
                // TODO: the value of FOLD_OVER_PENALTY need to be tested out
                edges[e_collapsed].error *= Internal::FOLD_OVER_PENALTY;
                heap.push(e_collapsed);
                continue;
            }

            // global index of kept edge in two deleted faces
            const vec2i e_kept{face2edge[ff[0]][v_del_in_ff[0]], face2edge[ff[1]][v_del_in_ff[1]]};

            // mark 2 faces, 1 vertex, and 1 edge as deleted;
            // 2 more edges will be marked in the loop later
            face_deleted[ff[0]] = true;
            face_deleted[ff[1]] = true;
            vertex_deleted[v_del] = true;
            edge_deleted[e_collapsed] = true;
            // update vertex position and quadric
            std::copy(edge.center.begin(), edge.center.end(), out_vertices[v_kept].begin());
            quadrics[v_kept] = edge.q;

            fve = fve_queue_v_del.front();
            Internal::Edge* tgt_edge = &edges[e_kept[0]];
            // first face to process: deleted face 0
            out_indices[f][3-v-e] = v_kept;
            edge_deleted[face2edge[f][e]] = true;
            face2edge[f][e] = e_kept[0];
            const idx ff0_in_edge = fi_in_edge(*tgt_edge, ff[0]);
            tgt_edge->faces[ff0_in_edge] = f;
            tgt_edge->idx_in_face[ff0_in_edge] = e;
            // every face centered around deleted vertex
            for (fve_queue_v_del.pop(); !fve_queue_v_del.empty(); fve_queue_v_del.pop()) {
                fve = fve_queue_v_del.front();
                out_indices[f][3-v-e] = v_kept;
                tgt_edge = &edges[face2edge[f][e]];
                Internal::update_edge(*tgt_edge, quadrics[out_indices[f][v]], edge.q,
                                      out_indices[f][v], v_kept);
            }
            // the deleted face 1
            edge_deleted[face2edge[f][v]] = true;
            face2edge[f][v] = e_kept[1];
            tgt_edge = &edges[e_kept[1]];
            const idx ff1_in_edge = fi_in_edge(*tgt_edge, ff[1]);
            tgt_edge->faces[ff1_in_edge] = f;
            tgt_edge->idx_in_face[ff1_in_edge] = v;

            // every face centered around the kept vertex
            for (; !fve_queue_v_kept.empty(); fve_queue_v_kept.pop()) {
                fve = fve_queue_v_kept.front();
                tgt_edge = &edges[face2edge[f][e]];
                Internal::update_edge(*tgt_edge, quadrics[out_indices[f][v]], edge.q);
            }
            tgt_edge = &edges[face2edge[f][v]];
            Internal::update_edge(*tgt_edge, quadrics[out_indices[f][e]], edge.q);
        }
        else {
            // TODO: handle boundary cases
            assert(false);
        }

        --nv;
    }

    // edges are pointless from this point on, but need to fix vertices and indices
    Internal::compact_data(out_vertices, out_indices, vertex_deleted, face_deleted);
    return {out_vertices, out_indices};
}

} // namespace MeshSimpl

#endif // MESH_SIMPL_LIBRARY_H
