//
// Created by nickl on 1/8/19.
//

#include "pre_proc.hpp"
#include "edge.hpp"
#include "util.hpp"
#include <limits>
#include <set>
#include <sstream>

namespace MeshSimpl {
namespace Internal {

void weight_quadric(Quadric& quadric, double face_area, WEIGHTING strategy) {
    switch (strategy) {
    case BY_AREA:
        quadric *= face_area;
        break;
    case BY_AREA_INV:
        if (face_area <= std::numeric_limits<double>::epsilon())
            quadric *= 0.0;
        else
            quadric *= 1 / face_area;
        break;
    case UNIFORM:
    default:
        if (face_area <= std::numeric_limits<double>::epsilon())
            quadric *= 0.0;
        break;
    }
}

Q compute_quadrics(const V& vertices, Internal::Connectivity& conn,
                   const SimplifyOptions& options) {
    // quadrics are initialized with all zeros
    Q quadrics(vertices.size());

    for (const auto& face : conn.indices) {
        // calculate the plane of this face (n and d: n'v+d=0 defines the plane)
        const vec3d edge01 = vertices[face[1]] - vertices[face[0]];
        const vec3d edge02 = vertices[face[2]] - vertices[face[0]];
        vec3d normal = cross(edge01, edge02);
        // |normal| = area, used for normalization and weighting quadrics
        const double area = magnitude(normal);
        normal /= area;
        // d = -n*v0
        const double d = -dot(normal, vertices[face[0]]);

        // calculate quadric Q = (A, b, c) = (nn', dn, d*d)
        Quadric q = make_quadric(normal, d);
        weight_quadric(q, area, options.weighting);

        for (const auto& v : face)
            quadrics[v] += q;
    }

    if (!options.fix_boundary) {
        for (idx f = 0; f < conn.indices.size(); ++f) {
            order k;
            for (k = 0; k < 3; ++k)
                if (conn.edge_of_face(f, k).on_boundary())
                    break;

            if (k == 3)
                continue;

            const auto& v = conn.indices[f];
            const std::array<vec3d, 3> e{vertices[v[2]] - vertices[v[1]],
                                         vertices[v[0]] - vertices[v[2]],
                                         vertices[v[1]] - vertices[v[0]]};
            const vec3d n_face = cross(e[0], e[1]);

            for (k = 0; k < 3; ++k) {
                if (!conn.edge_of_face(f, k).on_boundary())
                    continue;

                vec3d normal = cross(n_face, e[k]);
                const double normal_mag = magnitude(normal);
                normal /= normal_mag;
                const double d = -dot(normal, vertices[v[(k + 1) % 3]]);
                Quadric q = make_quadric(normal, d);
                q *= options.border_constraint;
                weight_quadric(q, magnitude(n_face), options.weighting);

                quadrics[v[(k + 1) % 3]] += q;
                quadrics[v[(k + 2) % 3]] += q;
            }
        }
    }

    return quadrics;
}

bool edge_topo_correctness(const Connectivity& conn) {
    for (idx f = 0; f < conn.indices.size(); ++f) {
        const auto& f2e = conn.face2edge[f];
        for (order i = 0; i < 3; ++i) {
            auto vv = conn.edges[f2e[i]].vertices;
            assert(vv[0] < vv[1]);
            idx v_smaller = conn.indices[f][(i + 1) % 3];
            idx v_larger = conn.indices[f][(i + 2) % 3];
            if (v_smaller > v_larger)
                std::swap(v_smaller, v_larger);

            // edge vertices should match face corners
            assert(vv[0] == v_smaller && vv[1] == v_larger);
        }
    }

    return true;
}

void construct_edges(size_t vertex_cnt, Internal::Connectivity& conn) {
    std::set<Edge> edge_set;

    // insert all edges into edge_set and find out if it is on boundary
    for (idx f = 0; f < conn.indices.size(); ++f) {
        const auto& face = conn.indices[f];
        for (order k = 0; k < 3; ++k) {
            // construct edge (v[i], v[j]);
            // edge local index will be k (= that of the 3rd vertex)
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
            edge.ord_in_faces[0] = k;
            edge.ord_in_faces[1] = Edge::INVALID;
            edge.boundary_v = Edge::BOTH;
            std::pair<std::set<Edge>::iterator, bool> it_and_inserted;
            it_and_inserted = edge_set.emplace(edge);
            auto it = it_and_inserted.first;

            if (!it_and_inserted.second) {
                if (it->boundary_v == Edge::NONE) {
                    std::ostringstream ss;
                    ss << "ERROR::NON_MANIFOLD_EDGE: face #" << (f + 1) << ", #"
                       << (it->faces[0] + 1) << ", and #" << (it->faces[1] + 1)
                       << " have an edge in common";
                    throw std::runtime_error(ss.str());
                }
                // edge was not inserted because it is already there
                // modifying through immutable iterator only if do not affect order
                *const_cast<Edge::BOUNDARY_V*>(&it->boundary_v) = Edge::NONE;
                *const_cast<idx*>(&it->faces[1]) = f;
                *const_cast<order*>(&it->ord_in_faces[1]) = k;
            }
        }
    }

    // convert edges from set to vector
    conn.edges = E(edge_set.begin(), edge_set.end());
    conn.face2edge.resize(conn.indices.size());
    std::vector<bool> vertex_on_boundary(vertex_cnt, false);

    for (idx i = 0; i < conn.edges.size(); ++i) {
        const auto& edge = conn.edges[i];
        // identify boundary vertices
        if (edge.boundary_v == Edge::BOTH) {
            vertex_on_boundary[edge.vertices[0]] = true;
            vertex_on_boundary[edge.vertices[1]] = true;
        }

        // populate face2edge references
        conn.face2edge[edge.faces[0]][edge.ord_in_faces[0]] = i;
        if (edge.boundary_v != Edge::BOTH)
            conn.face2edge[edge.faces[1]][edge.ord_in_faces[1]] = i;
    }

    // non-boundary edges may have one vertex on boundary, find them in this loop
    for (auto& edge : conn.edges) {
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

    assert(edge_topo_correctness(conn));
}

} // namespace Internal
} // namespace MeshSimpl
