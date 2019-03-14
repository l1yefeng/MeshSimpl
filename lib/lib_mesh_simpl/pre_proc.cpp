//
// Created by nickl on 1/8/19.
//

#include "pre_proc.h"
#include "util.h"
#include <set>
#include <sstream>

namespace MeshSimpl {

namespace Internal {

void compute_quadrics_per_face(const V& vertices, const std::vector<idx>& face,
                               const bool weight_by_area, Q& quadrics) {
    // calculate the plane of this face (n and d: n'v+d=0 defines the plane)
    const vec3d edge1 = vertices[face[1]] - vertices[face[0]];
    const vec3d edge2 = vertices[face[2]] - vertices[face[0]];
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

Q compute_quadrics(const V& vertices, const F& indices, bool weight_by_area) {
    // quadrics are initialized with all zeros
    Q quadrics(vertices.size());
    for (const auto& face : indices)
        compute_quadrics_per_face(vertices, face, weight_by_area, quadrics);
    return quadrics;
}

Q recompute_quadrics(const V& vertices, const F& indices, const std::vector<bool>& deleted_face,
                     bool weight_by_area) {
    // quadrics are initialized with all zeros
    Q quadrics(vertices.size());
    for (idx f = 0; f < indices.size(); ++f)
        if (!deleted_face[f])
            compute_quadrics_per_face(vertices, indices[f], weight_by_area, quadrics);
    return quadrics;
}

bool edge_topo_correctness(const E& edges, const F2E& face2edge, const F& indices) {
    for (idx f = 0; f < indices.size(); ++f) {
        const auto& f2e = face2edge[f];
        for (order i = 0; i < 3; ++i) {
            auto vv = edges[f2e[i]].vertices;
            if (vv[0] >= vv[1])
                return false;
            idx v_smaller = indices[f][(i + 1) % 3];
            idx v_larger = indices[f][(i + 2) % 3];
            if (v_smaller > v_larger)
                std::swap(v_smaller, v_larger);

            // edge vertices should match face corners
            if (vv[0] != v_smaller || vv[1] != v_larger)
                return false;
        }
    }

    return true;
}

std::pair<E, F2E> construct_edges(const F& indices, const size_t vertex_cnt) {
    const auto edge_cmp = [](const Edge& a, const Edge& b) -> bool {
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
        for (order k = 0; k < 3; ++k) {
            // construct edge (v[i], v[j]); edge local index will be k (= that of the 3rd vertex)
            const order i = (k + 1) % 3;
            const order j = (k + 2) % 3;
            idx v0 = face[i];
            idx v1 = face[j];
            if (v0 > v1)
                std::swap(v0, v1);
            Edge edge{};
            edge.vertices[0] = v0;
            edge.vertices[1] = v1;
            edge.faces[0] = f;
            edge.idx_in_face[0] = k;
            edge.boundary_v = Edge::BOTH;
            auto it_and_inserted = edge_set.emplace(edge);
            auto it = it_and_inserted.first;
            if (!it_and_inserted.second) {
                if (it->boundary_v == Edge::NONE) {
                    std::ostringstream ss;
                    ss << "ERROR::NON_MANIFOLD_EDGE: please check face #" << f << ", #"
                       << it->faces[0] << ", and #" << it->faces[1];
                    throw std::runtime_error(ss.str());
                }
                // edge was not inserted because it is already there
                // modifying through immutable iterator only if do not affect order
                *const_cast<Edge::BOUNDARY_V*>(&it->boundary_v) = Edge::NONE;
                *const_cast<idx*>(&it->faces[1]) = f;
                *const_cast<order*>(&it->idx_in_face[1]) = k;
            }
        }
    }

    // convert edges from set to vector
    E edges(edge_set.begin(), edge_set.end());
    std::vector<bool> vertex_on_boundary(vertex_cnt, false);
    F2E face2edge(indices.size());

    for (idx i = 0; i < edges.size(); ++i) {
        const auto& edge = edges[i];
        // identify boundary vertices
        if (edge.boundary_v == Edge::BOTH) {
            vertex_on_boundary[edge.vertices[0]] = true;
            vertex_on_boundary[edge.vertices[1]] = true;
        }

        // populate face2edge references
        face2edge[edge.faces[0]][edge.idx_in_face[0]] = i;
        if (edge.boundary_v != Edge::BOTH)
            face2edge[edge.faces[1]][edge.idx_in_face[1]] = i;
    }

    // non-boundary edges may have one vertex on boundary, find them in this loop
    for (auto& edge : edges) {
        if (edge.boundary_v == Edge::BOTH)
            continue;
        if (vertex_on_boundary[edge.vertices[0]] && vertex_on_boundary[edge.vertices[1]])
            edge.boundary_v = Edge::BOTH;
        else if (vertex_on_boundary[edge.vertices[0]])
            edge.boundary_v = Edge::V0;
        else if (vertex_on_boundary[edge.vertices[1]])
            edge.boundary_v = Edge::V1;
        // else totally within the boundary
    }

    assert(edge_topo_correctness(edges, face2edge, indices));

    return {edges, face2edge};
}

} // namespace Internal

} // namespace MeshSimpl
