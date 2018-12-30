#ifndef MESH_SIMPL_LIBRARY_H
#define MESH_SIMPL_LIBRARY_H

#include "util.h"
#include <vector>
#include <array>
#include <set>
#include <ostream>
#include <queue>
#include <algorithm>

namespace MeshSimpl
{

namespace Internal
{

enum BOUNDARY_V { NONE, BOTH, V0, V1 };

// Edge defines the struct of an edge
struct Edge
{
    std::array<idx, 2> vertices;    // index of two end vertices, unique, v0 < v1
    std::array<idx, 2> faces;       // index of two incident faces
    std::array<idx, 2> idx_in_face; // index of this (0, 1, 2) in faces
    BOUNDARY_V boundary_v;          // vertices on boundary
    vec3d center;                   // where this edge collapse into
    double error;                   // quadric error value, undefined if on boundary
    Quadric q;                      // sum of quadrics of two vertices
    bool dirty;                     // true if this edge has wrong position in heap
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
std::pair<std::vector<Edge>, std::vector<vec3u>>
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
            idx j = (i+1)%3;
            idx v0 = face[i], v1 = face[j];
            if (v0 > v1)
                std::swap(v0, v1);
            auto it_and_inserted = edge_set.insert({{v0, v1}, {f}, {i}, BOUNDARY_V::BOTH});
            auto it = it_and_inserted.first;
            if (!it_and_inserted.second) {
                // edge was not inserted because it is already there
                // modifying through immutable iterator only if do not affect order
                *const_cast<BOUNDARY_V*>(&it->boundary_v) = BOUNDARY_V::NONE;
                *const_cast<idx*>(&it->faces[1]) = f;
                *const_cast<idx*>(&it->idx_in_face[1]) = (i+2)%3;
            }
        }
    }

    // convert edges from set to vector
    std::vector<Edge> edges(edge_set.begin(), edge_set.end());
    std::vector<bool> vertex_on_boundary(vertex_cnt, false);
    std::vector<vec3u> face2edge(indices.size());

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
    size_t i = 0, j = indices.size()-1;
    while (true) {
        while (!face_deleted[i] && i <= j)
            ++i;
        while (face_deleted[j] && i < j)
            --j;
        if (i >= j)
            break;
        std::swap(indices[i++], indices[j--]);
    }
    indices.resize(i);

    // create mapping from v to f
    for (i = 0; i < vertices.size(); ++i)
        for (j = 0; j < 3; ++j)
            vertex2face[indices[i][j]][j].push_back(static_cast<idx>(i));

    // get rid of deleted vertices and keep the mapping valid
    i = 0;
    j = vertices.size()-1;
    while (true) {
        while (!vertex_deleted[i] && i <= j)
            ++i;
        while (vertex_deleted[j] && i < j)
            --j;
        if (i >= j)
            break;
        std::swap(vertices[i], vertices[j]);
        std::swap(vertex2face[i], vertex2face[j]);
        ++i;
        --j;
    }
    vertices.resize(i);
    vertex2face.resize(i);

    for (i = 0; i < vertex2face.size(); ++i)
        for (j = 0; j < 3; ++j)
            for (const idx f : vertex2face[i][j])
                indices[f][j] = static_cast<idx>(i);
}

} // namespace MeshSimpl::Internal

// Mesh simplification main method. Simplify given mesh until remaining number of vertices/faces
// is (1-strength) of the original. Returns output vertices and indices as in inputs.
// TODO: Ignoring the 4th and the following values (if exist) in vertices.
std::pair<V, F> simplify(const V& vertices, const F& indices, float strength)
{
    const auto out_nv = static_cast<size_t>(std::lround((1-strength)*vertices.size()));
    auto out_vertices = vertices;
    auto out_indices = indices;

    auto quadrics = Internal::compute_quadrics(vertices, indices, true);

    auto edge_topo = Internal::edge_topology(indices, vertices.size());
    auto& edges = edge_topo.first;
    auto& face2edge = edge_topo.second;

    Internal::init_edge_errors(vertices, quadrics, edges);

    // create priority queue on quadric error
    auto heap = Internal::build_min_heap(edges);
    size_t nv = vertices.size();
    std::vector<bool> vertex_deleted(vertices.size(), false);
    std::vector<bool> face_deleted(indices.size(), false);
    std::vector<bool> edge_deleted(edges.size(), false);

    // returns the position of vertex in face. a convenient method but make sure v EXISTS in face
    const auto vi_in_face = [&](const idx v, const idx f)->idx {
        const auto& vec = out_indices[f];
        return static_cast<idx>(std::distance(vec.begin(), std::find(vec.begin(), vec.end(), v)));
    };

    // returns f,v,e of next neighbor face centered around v_center;
    // indexes v and e are local in face (0,1,2) in both input and output.
    const auto iter_next = [&](const vec3u& fve)->vec3u {
        const idx f = fve[0], v = fve[1], e = fve[2];
        const auto& edge = edges[face2edge[f][v]];
        const idx f_idx_to_edge = edge.faces[0] == f ? 0 : 1;
        const idx of = edge.faces[1-f_idx_to_edge];
        const idx oe = edge.idx_in_face[1-f_idx_to_edge];
        const idx ov_global = out_indices[f][e];
        const idx ov = vi_in_face(of, ov_global);
        return {of, ov, oe};
    };

    while (!heap.empty() && nv > out_nv) {
        // collapse an edge to remove 1 vertex and 2 faces in each iteration
        const idx e_collapsed = heap.top();
        heap.pop();

        if (edge_deleted[e_collapsed])
            continue;

        const auto& edge = edges[e_collapsed];
        if (edge.dirty) {
            // this edge was modified and is having an error no less than before
            edges[e_collapsed].dirty = true;
            heap.push(e_collapsed);
            continue;
        }

        // some aliases
        const auto& vv = edge.vertices;
        const auto& ff = edge.faces;

        if (edge.boundary_v == Internal::BOUNDARY_V::NONE) {
            // the choice of deletion does not matter
            const idx v_kept = vv[0];
            const idx v_del = vv[1];
            // mark as deleted
            face_deleted[ff[0]] = true;
            face_deleted[ff[1]] = true;
            vertex_deleted[v_del] = true;
            // update vertex position and quadric
            std::copy(edge.center.begin(), edge.center.end(), out_vertices[v_kept].begin());
            quadrics[v_kept] = edge.q;

            // starting from ff[0], iterate through faces centered around v_del and modify
            std::array<idx, 2> v_kept_in_ff{vi_in_face(ff[0], v_kept), vi_in_face(ff[1], v_kept)};
            std::array<idx, 2> v_del_in_ff{vi_in_face(ff[0], v_del), vi_in_face(ff[1], v_del)};
            vec3u fve = iter_next({ff[0], v_kept_in_ff[0], edge.idx_in_face[0]});
            vec3u fve_next;
            std::array<idx, 2> e_kept{};
            for (idx i = 0; i < 2; ++i)
                e_kept[i] = face2edge[ff[i]][v_del_in_ff[i]];

            bool first_iter = true;
            while (true) {
                const idx f = fve[0], v = fve[1], e = fve[2];
                fve_next = iter_next(fve);
                // v_del got replaced by v_kept
                out_indices[f][vi_in_face(v_del, f)] = v_kept;
                // deal with the first edge
                if (first_iter) {
                    edge_deleted[face2edge[f][e]] = true;
                    face2edge[f][e] = e_kept[0];

                    // update faces and idx_in_faces of the edge in between
                    auto& edge_kept = edges[e_kept[0]];
                    const idx fi_in_edge = edge_kept.faces[0] == ff[0] ? 0 : 1;
                    edge_kept.faces[fi_in_edge] = f;
                    edge_kept.idx_in_face[fi_in_edge] = e;
                }
                else {
                    auto& edge_dirty = edges[face2edge[f][e]];
                    Internal::update_edge(edge_dirty, quadrics[out_indices[f][v]], edge.q,
                                          out_indices[f][v], v_kept);
                }
                // iteration will end because next face is deleted
                if (fve_next[0] == ff[1]) {
                    edge_deleted[face2edge[f][v]] = true;
                    face2edge[f][v] = e_kept[1];

                    // update faces and idx_in_faces of the edge in between
                    auto& edge_kept = edges[e_kept[1]];
                    const idx fi_in_edge = edge_kept.faces[0] == ff[1] ? 0 : 1;
                    edge_kept.faces[fi_in_edge] = f;
                    edge_kept.idx_in_face[fi_in_edge] = v;
                    break;
                }

                fve = fve_next;
                first_iter = false;
            }

            // starting from ff[1], iterate through faces centered around v_kept and update edges
            fve = iter_next({fve_next[0], v_del_in_ff[1], edge.idx_in_face[1]});
            for (;; fve = iter_next(fve)) {
                const idx f = fve[0], v = fve[1], e = fve[2];
                auto& edge_dirty = edges[face2edge[f][e]];
                Internal::update_edge(edge_dirty, quadrics[out_indices[f][v]], edge.q);
                if (face2edge[f][e] == e_kept[0])
                    break;
            }
        }
        else {
            // TODO: handle boundary cases
        }

        --nv;
    }

    // edges are pointless from this point on, but need to fix vertices and indices
    Internal::compact_data(out_vertices, out_indices, vertex_deleted, face_deleted);
    return {out_vertices, out_indices};
}

} // namespace MeshSimpl

#endif // MESH_SIMPL_LIBRARY_H
