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

bool is_valid_edge_target(const Internal::Edge& edge,
                          const std::vector<bool>& deleted_face,
                          const std::vector<bool>& deleted_vertex, bool fix_boundary) {
    assert(edge.vertices[0] != edge.vertices[1]);
    assert(!deleted_vertex[edge.vertices[0]] && !deleted_vertex[edge.vertices[1]]);
    if (fix_boundary) {
        assert(edge.boundary_v != Internal::Edge::BOTH);
    }

    if (edge.on_boundary()) {
        assert(!deleted_face[edge.faces[0]]);
        assert(edge.boundary_v == Internal::Edge::BOTH);
    } else {
        assert(edge.faces[0] != edge.faces[1]);
        assert(!deleted_face[edge.faces[0]] && !deleted_face[edge.faces[1]]);
    }
    return true;
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

    std::vector<char> boundary_flags;
    if (!options.fix_boundary)
        boundary_flags = std::vector<char>(indices.size(), 0x07);

    // [1] find out information of edges (endpoints, incident faces) and face2edge
    auto edge_topo = Internal::construct_edges(indices, NV, boundary_flags);
    std::vector<Internal::Edge>& edges = edge_topo.first;
    std::vector<vec3i>& face2edge = edge_topo.second;

    // [2] compute quadrics of vertices
    auto quadrics =
        Internal::compute_quadrics(vertices, indices, boundary_flags, options.weighting);

    if (!boundary_flags.empty()) {
        // throw away
        boundary_flags = std::vector<char>();
    }

    // [3] assigning edge errors using quadrics
    Internal::compute_errors(vertices, quadrics, edges, options.fix_boundary);

    // [4] create priority queue on quadric error
    Internal::QEMHeap heap(edges, !options.fix_boundary);

    std::vector<bool> deleted_face(indices.size(), false);
    std::vector<bool> deleted_vertex(NV, true);
    for (const auto& face : indices)
        for (idx v : face)
            deleted_vertex[v] = false;

    size_t nv = nv_to_decimate;
    idx prev_target = static_cast<idx>(edges.size());
    while (!heap.empty() && nv > 0) {
        // target the least-error edge, if it is what we saw last iteration,
        // it means loop should stop because all remaining edges have been penalized
        const idx target = heap.top();
        if (prev_target == target)
            break;
        prev_target = target;

        const auto& edge = edges[target];

        assert(is_valid_edge_target(edge, deleted_face, deleted_vertex,
                                    options.fix_boundary));

        // [5] collapse the least-error edge until mesh is simplified enough
        bool collapsed = edge_collapse(out_vertices, out_indices, edges, face2edge,
                                       quadrics, heap, target, options.fix_boundary);
        if (!collapsed)
            continue;

        // mark adjacent faces deleted
        deleted_face[edge.faces[0]] = true;
        if (edge.boundary_v != Internal::Edge::BOTH)
            deleted_face[edge.faces[1]] = true;
        else
            assert(edge.on_boundary());

        // mark one (chosen) endpoint deleted
        deleted_vertex[edge.vertices[edge.v_del_order()]] = true;

        // of course there might be edges deleted,
        // they must have been removed from heap during `edge_collapse`

        --nv;
    }

    // edges are pointless from this point on, but need to fix vertices and indices
    Internal::compact_data(deleted_vertex, deleted_face, out_vertices, out_indices);

    return {out_vertices, out_indices};
}

} // namespace MeshSimpl
