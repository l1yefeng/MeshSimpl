//
// Created by nickl on 1/8/19.
//

#include "simplify.hpp"
#include "ecol.hpp"
#include "marker.hpp"
#include "post_proc.hpp"
#include "pre_proc.hpp"
#include <limits>

namespace MeshSimpl {

using namespace Internal;

void options_validation(const SimplifyOptions& options) {
    if (options.strength > 1)
        throw std::invalid_argument("ERROR::INVALID_OPTION: strength > 1");
    if (options.strength < 0)
        throw std::invalid_argument("ERROR::INVALID_OPTION: strength < 0");
}

bool is_valid_edge_target(const Edge& edge, const Marker& marker,
                          bool fix_boundary) {
    assert(edge.vertices[0] != edge.vertices[1]);
    assert(marker.exist_v(edge.vertices[0]) &&
           marker.exist_v(edge.vertices[1]));
    if (fix_boundary) {
        assert(edge.boundary_v != Edge::BOTH);
    }

    if (edge.on_boundary()) {
        assert(marker.exist_f(edge.faces[0]));
        assert(edge.boundary_v == Edge::BOTH);
    } else {
        assert(edge.faces[0] != edge.faces[1]);
        assert(marker.exist_f(edge.faces[0]) && marker.exist_f(edge.faces[1]));
    }
    return true;
}

bool simplify(TriMesh& mesh, const SimplifyOptions& options) {
    options_validation(options);

    const size_t NV = mesh.vertices.size();
    const size_t nv_to_decimate =
        NV -
        std::max(static_cast<int>(std::lround((1 - options.strength) * NV)), 3);
    auto& vertices = mesh.vertices;
    auto& indices = mesh.indices;

    if (nv_to_decimate == 0)
        return true;

    // keeps a reference to indices and init face2edge and edges
    Connectivity conn{indices, {}, {}};

    // marker keeps track of face/vertex state: exist or deleted
    Marker marker(NV, indices.size());
    marker.mark_unref_v(indices);

    // [1] find out information of edges (endpoints, incident faces) and face2edge
    construct_edges(vertices, conn);

    // [2] compute quadrics of vertices
    Q quadrics;
    compute_quadrics(vertices, conn, quadrics, options);

    // [3] assigning edge errors using quadrics
    compute_errors(vertices, quadrics, conn.edges, options.fix_boundary);

    // [4] create priority queue on quadric error
    QEMHeap heap(conn.edges, !options.fix_boundary);

    size_t nv = nv_to_decimate;
    while (!heap.empty() && nv > 0) {
        // target the least-error edge, if it is what we saw last iteration,
        // it means loop should stop because all remaining edges have been penalized
        const idx target = heap.top();
        if (conn.edges[target].error >= std::numeric_limits<double>::max())
            break;

        const auto& edge = conn.edges[target];

        assert(is_valid_edge_target(edge, marker, options.fix_boundary));

        // [5] collapse the least-error edge until mesh is simplified enough
        bool collapsed =
            edge_collapse(vertices, conn, quadrics, heap, target, options);
        if (!collapsed)
            continue;

        // mark adjacent faces deleted
        marker.mark_f(edge.faces[0]);
        if (edge.boundary_v != Edge::BOTH)
            marker.mark_f(edge.faces[1]);
        else
            assert(edge.on_boundary());

        // mark one (chosen) endpoint deleted
        marker.mark_v(edge.vertices[edge.v_del_order()]);

        // of course there might be edges deleted,
        // they must have been removed from heap during `edge_collapse`

        --nv;
    }

    // edges are pointless from this point on, but need to fix vertices and indices
    compact_data(vertices, indices, marker);

    return true;
}

} // namespace MeshSimpl
