//
// Created by nickl on 1/8/19.
//

#include "pre_proc.h"
#include "util.h"
#include <set>

namespace MeshSimpl
{
namespace Internal
{

void compute_quadrics_per_face(const V& vertices, const std::vector<idx>& face,
                               const bool weight_by_area, std::vector<Quadric>& quadrics)
{
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

std::vector<Quadric> compute_quadrics(const V& vertices, const F& indices, bool weight_by_area)
{
    // quadrics are initialized with all zeros
    std::vector<Quadric> quadrics(vertices.size());
    for (const auto& face : indices)
        compute_quadrics_per_face(vertices, face, weight_by_area, quadrics);
    return quadrics;
}

std::vector<Quadric>
compute_quadrics(const V& vertices, const F& indices, const std::vector<bool>& deleted_face,
                 bool weight_by_area)
{
    // quadrics are initialized with all zeros
    std::vector<Quadric> quadrics(vertices.size());
    for (idx f = 0; f < indices.size(); ++f)
        if (!deleted_face[f])
            compute_quadrics_per_face(vertices, indices[f], weight_by_area, quadrics);
    return quadrics;
}

std::pair<std::vector<Edge>, std::vector<vec3i>>
construct_edges(const F& indices, const size_t vertex_cnt)
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

}
}