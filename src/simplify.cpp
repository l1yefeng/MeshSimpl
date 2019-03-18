//
// Created by nickl on 1/8/19.
//

#include "simplify.hpp"
#include "ecol.hpp"
#include "post_proc.hpp"
#include "pre_proc.hpp"

namespace MeshSimpl {

void options_validation(const SimplifyOptions& options) {
    if (options.strength > 1)
        throw std::invalid_argument("ERROR::INVALID_OPTION: strength > 1");
    if (options.strength < 0)
        throw std::invalid_argument("ERROR::INVALID_OPTION: strength < 0");
}

std::pair<V, F> simplify(const V& vertices, const F& indices,
                         const SimplifyOptions& options) {
    options_validation(options);

    const size_t NV = vertices.size();
    const size_t nv_to_decimate =
        NV - std::max(static_cast<int>(std::lround((1 - options.strength) * NV)), 4);
    auto out_vertices = vertices;
    auto out_indices = indices;

    if (nv_to_decimate == 0)
        return {out_vertices, out_indices};

    Internal::ConstraintPlane constraint_plane;
    if (options.fix_boundary) {
        constraint_plane.enabled = false;
    } else {
        constraint_plane.enabled = true;
        constraint_plane.on_boundary.resize(indices.size(), true);
        constraint_plane.boundary_e_order.resize(indices.size());
    }

    // [1] find out information of edges (endpoints, incident faces) and face2edge
    auto edge_topo = Internal::construct_edges(indices, NV, constraint_plane);
    std::vector<Internal::Edge>& edges = edge_topo.first;
    std::vector<vec3i>& face2edge = edge_topo.second;

    // [2] compute quadrics of vertices
    auto quadrics = Internal::compute_quadrics(vertices, indices, constraint_plane, options.weighting);

    if (constraint_plane.enabled) {
        // throw away
        constraint_plane.on_boundary = std::vector<bool>();
        constraint_plane.boundary_e_order = std::vector<order>();
    }

    // [3] assigning edge errors using quadrics
    Internal::compute_errors(vertices, quadrics, edges, options.fix_boundary);

    // [4] create priority queue on quadric error
    Internal::QEMHeap heap(edges, !options.fix_boundary);

    std::vector<bool> deleted_face(indices.size(), false);
    std::vector<bool> deleted_vertex(NV, true);
    for (const auto& face : indices) {
        for (idx v : face) {
            deleted_vertex[v] = false;
        }
    }

    size_t nv = nv_to_decimate;
    idx prev_target = static_cast<idx>(edges.size());
    while (!heap.empty() && nv > 0) {
        // collapse an edge to remove 1 vertex and 2 faces in each iteration
        const idx target = heap.top();
        if (prev_target == target)
            break;
        prev_target = target;

        const auto& edge = edges[target];

        assert(edge.vertices[0] != edge.vertices[1]);
        assert(!deleted_face[edge.faces[0]]);
        assert(!deleted_vertex[edge.vertices[0]] && !deleted_vertex[edge.vertices[1]]);

        if (options.fix_boundary) {
            assert(edge.faces[0] != edge.faces[1]);
            assert(edge.boundary_v != Internal::Edge::BOTH);
            assert(!deleted_face[edge.faces[1]]);
        }

        // [5] collapse the least-error edge until mesh is simplified enough
        bool collapsed = edge_collapse(out_vertices, out_indices, edges, face2edge,
                                       quadrics, heap, target, options.fix_boundary);
        if (!collapsed)
            continue;

        deleted_face[edge.faces[0]] = true;
        if (edge.boundary_v != Internal::Edge::BOTH)
            deleted_face[edge.faces[1]] = true;

        deleted_vertex[edge.vertices[Internal::choose_v_del(edge)]] = true;
        --nv;
    }

    // edges are pointless from this point on, but need to fix vertices and indices
    Internal::compact_data(deleted_vertex, deleted_face, out_vertices, out_indices);

    return {out_vertices, out_indices};
}

} // namespace MeshSimpl
