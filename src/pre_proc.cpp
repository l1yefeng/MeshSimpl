//
// Created by nickl on 1/8/19.
//

#include "pre_proc.hpp"
#include "edge.hpp"
#include "util.hpp"
#include <algorithm>
#include <limits>

namespace MeshSimpl {
namespace Internal {

void weight_quadric(Quadric& quadric, double face_area, WEIGHTING strategy) {
    assert(face_area != 0);
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

void compute_quadrics(const V& vertices, Internal::Connectivity& conn, Q& quadrics,
                      const std::vector<bool>& deleted_face,
                      const SimplifyOptions& options) {
    // quadrics are initialized with all zeros
    quadrics.resize(vertices.size());

    for (idx f = 0; f < conn.indices.size(); ++f) {
        if (deleted_face[f])
            continue;

        const auto& face = conn.indices[f];
        // calculate the plane of this face (n and d: n'v+d=0 defines the plane)
        const vec3d edge01 = vertices[face[1]] - vertices[face[0]];
        const vec3d edge02 = vertices[face[2]] - vertices[face[0]];
        vec3d normal = cross(edge01, edge02);
        // |normal| = area, used for normalization and weighting quadrics
        const double area = magnitude(normal);
        if (area != 0)
            normal /= area;
        else
            continue;

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
            if (deleted_face[f])
                continue;

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
                if (normal_mag != 0)
                    normal /= normal_mag;
                else
                    continue;

                const double d = -dot(normal, vertices[v[(k + 1) % 3]]);
                Quadric q = make_quadric(normal, d);
                q *= options.border_constraint;
                weight_quadric(q, magnitude(n_face), options.weighting);

                quadrics[v[(k + 1) % 3]] += q;
                quadrics[v[(k + 2) % 3]] += q;
            }
        }
    }
}

bool edge_topo_correctness(const Connectivity& conn,
                           const std::forward_list<idx>& invalid_faces) {
    for (idx f = 0; f < conn.indices.size(); ++f) {
        if (flist_contains(invalid_faces, f))
            continue;

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

void ensure_manifold(std::map<std::pair<idx, idx>, std::pair<Edge, flist_t>>& edge_set,
                     const V& vertices, const F& indices,
                     std::forward_list<idx>& invalid_faces) {
    for (auto& elem : edge_set) {
        const auto& face_list = elem.second.second;
        // non-manifolds have non-empty face_list
        if (face_list.empty())
            continue;

        auto& edge = elem.second.first;

        // will keep faces with largest area valid, thus put area first to sort
        typedef std::tuple<double, idx, order> candidate_t;
        std::vector<candidate_t> candidates;
        candidates.reserve(4);
        for (const std::pair<idx, order>& fo : face_list) {
            const auto& face = indices[fo.first];
            // if this face is already invalid because of some previous non-manifold edge
            if (flist_contains(invalid_faces, fo.first))
                continue;

            candidates.emplace_back(std::make_tuple(
                tri_area(vertices[face[0]], vertices[face[1]], vertices[face[2]]),
                fo.first, fo.second));
        }

        if (candidates.empty()) {
            // the edge have been non-manifold, but some other edge removed some face(s)
            // and now it is not any more
            continue;
        }

        assert(edge.ord_in_faces[0] != Edge::INVALID &&
               edge.ord_in_faces[1] != Edge::INVALID);
        for (order i : {0, 1}) {
            const auto& face = indices[edge.faces[i]];
            candidates.emplace_back(std::make_tuple(
                tri_area(vertices[face[0]], vertices[face[1]], vertices[face[2]]),
                edge.faces[i], edge.ord_in_faces[i]));
        }

        // update edge and drop the others
        std::sort(candidates.begin(), candidates.end(), std::greater<candidate_t>());
        for (order i : {0, 1}) {
            edge.faces[i] = std::get<1>(candidates[i]);
            edge.ord_in_faces[i] = std::get<2>(candidates[i]);
        }

        // about to drop some faces but they are referenced by some edges currently
        for (auto it = candidates.begin() + 2; it != candidates.end(); ++it) {
            idx f = std::get<1>(*it);
            order ord = std::get<2>(*it);
            invalid_faces.emplace_front(f);

            for (auto i : {1, 2}) {
                idx v0 = indices[f][ord];
                idx v1 = indices[f][(ord + i) % 3];
                if (v0 > v1)
                    std::swap(v0, v1);
                auto set_iter = edge_set.find(std::make_pair(v0, v1));
                assert(set_iter != edge_set.end());

                Edge& cur_edge = set_iter->second.first;
                auto& cur_flist = set_iter->second.second;

                // pop invalid faces from the front
                // most of the time it does nothing but when it does, filling
                // the vacant face would be safely operated (in else statement)
                while (!cur_flist.empty() &&
                       flist_contains(invalid_faces, cur_flist.front().first))
                    cur_flist.pop_front();

                // this invalid face is in the flist of `cur_edge`
                if (cur_edge.faces[0] != f && cur_edge.faces[1] != f) {
                    assert(std::find_if(cur_flist.begin(), cur_flist.end(),
                                        [&](const std::pair<idx, order>& x) -> bool {
                                            return x.first == f;
                                        }) != cur_flist.end());
                    continue;
                }

                const order f_in_edge = cur_edge.f_order(f);

                if (cur_flist.empty()) {
                    // normal edge becomes boundary edge
                    cur_edge.ord_in_faces[f_in_edge] = Edge::INVALID;
                    cur_edge.faces[f_in_edge] = 0; // does not really matter but anyway

                    if (cur_edge.ord_in_faces[1 - f_in_edge] == Edge::INVALID) {
                        // was boundary edge, now should disappear
                        assert(cur_edge.faces[1 - f_in_edge] == 0);
                    } else {
                        // was normal 2-manifold edge, now should be on boundary
                        assert(cur_edge.boundary_v == Edge::NEITHER);
                        cur_edge.boundary_v = Edge::BOTH;
                        if (f_in_edge == 0)
                            cur_edge.swap_faces();
                    }
                } else {
                    // a face in the flist should fill the vacant of cur_edge.faces
                    assert(cur_edge.boundary_v == Edge::NEITHER);
                    std::pair<idx, order>& replacement = cur_flist.front();
                    cur_edge.faces[f_in_edge] = replacement.first;
                    cur_edge.ord_in_faces[f_in_edge] = replacement.second;
                    cur_flist.pop_front();
                }
            }
        }
    }
}

void construct_edges(const V& vertices, Internal::Connectivity& conn,
                     std::forward_list<idx>& invalid_faces) {
    std::map<std::pair<idx, idx>, std::pair<Edge, flist_t>> edge_set;

    bool non_manifold_exist = false;

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
            auto it_and_inserted =
                edge_set.emplace(std::make_pair(v0, v1), std::make_pair(edge, flist_t{}));
            auto it = it_and_inserted.first;

            if (!it_and_inserted.second) {
                auto& curr_edge = it->second.first;
                if (curr_edge.boundary_v == Edge::NEITHER) {
                    // add to face list
                    it->second.second.emplace_front(std::make_pair(f, k));
                    non_manifold_exist = true;
                    continue;
                }
                curr_edge.boundary_v = Edge::NEITHER;
                curr_edge.faces[1] = f;
                curr_edge.ord_in_faces[1] = k;
            }
        }
    }

    if (non_manifold_exist) {
        // handle non-manifold edges, keep index of deleted faces in the list
        ensure_manifold(edge_set, vertices, conn.indices, invalid_faces);
    }

    // populate edges vector from map
    conn.edges.reserve(edge_set.size());
    if (!non_manifold_exist) {
        for (const auto& elem : edge_set)
            conn.edges.emplace_back(elem.second.first);
    } else {
        for (const auto& elem : edge_set)
            if (elem.second.first.ord_in_faces[0] != Edge::INVALID)
                conn.edges.emplace_back(elem.second.first);
    }

    conn.face2edge.resize(conn.indices.size());
    std::vector<bool> vertex_on_boundary(vertices.size(), false);

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

    assert(edge_topo_correctness(conn, invalid_faces));
}

} // namespace Internal
} // namespace MeshSimpl
