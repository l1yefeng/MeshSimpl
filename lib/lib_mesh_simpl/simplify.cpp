//
// Created by nickl on 1/8/19.
//

#include "simplify.h"
#include "ecol.h"
#include "post_proc.h"
#include "pre_proc.h"
#include <limits>

namespace MeshSimpl {

void options_validation(const SimplifyOptions& options) {
    if (options.strength > 1)
        throw std::invalid_argument("ERROR::INVALID_OPTION: strength > 1");
    if (options.strength < 0)
        throw std::invalid_argument("ERROR::INVALID_OPTION: strength < 0");
}

std::pair<V, F> simplify(const V& vertices, const F& indices, const SimplifyOptions& options) {
    options_validation(options);
    const size_t NV = vertices.size();
    size_t nv_to_decimate =
        NV - std::max(static_cast<size_t>(std::lround((1 - options.strength) * NV)),
                      Internal::MIN_NR_VERTICES);
    auto out_vertices = vertices;
    auto out_indices = indices;

    if (nv_to_decimate == 0)
        return {out_vertices, out_indices};

    auto quadrics = Internal::compute_quadrics(vertices, indices, options.weight_by_face);

    auto edge_topo = Internal::construct_edges(indices, NV);
    std::vector<Internal::Edge>& edges = edge_topo.first;
    std::vector<vec3i>& face2edge = edge_topo.second;

    Internal::compute_errors(vertices, quadrics, edges);

    // create priority queue on quadric error
    Internal::QEMHeap heap(edges);

    std::vector<bool> deleted_vertex(NV, false);
    std::vector<bool> deleted_face(indices.size(), false);

    // one run of a series of edge collapse that is supposed to decimate nv_decimate vertices and
    // returns true if it ends because this target is achieved instead of because of other reason
    while (!heap.empty() && nv_to_decimate > 0) {
        // collapse an edge to remove 1 vertex and 2 faces in each iteration
        const idx target = heap.top();

        const auto& edge = edges[target];
        if (edge.error == std::numeric_limits<double>::max())
            break;
        assert(edge.vertices[0] != edge.vertices[1]);
        assert(edge.faces[0] != edge.faces[1]);

        assert(edge.boundary_v != Internal::BOUNDARY_V::BOTH);
        if (edge_collapse(out_vertices, out_indices, edges, face2edge, quadrics, heap, target)) {
            for (const idx f : edge.faces)
                deleted_face[f] = true;
            deleted_vertex[edge.vertices[Internal::choose_v_del(edge)]] = true;
            --nv_to_decimate;
        }
    }

    // edges are pointless from this point on, but need to fix vertices and indices
    Internal::compact_data(deleted_vertex, deleted_face, out_vertices, out_indices);
    return {out_vertices, out_indices};
}

} // namespace MeshSimpl
