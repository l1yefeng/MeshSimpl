//
// Created by nickl on 1/8/19.
//

#include "simplify.hpp"
#include "ecol.hpp"
#include "post_proc.hpp"
#include "pre_proc.hpp"
#include <limits>

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

bool simplify(TriMesh& mesh, const SimplifyOptions& options) {
    options_validation(options);

    const size_t NV = mesh.vertices.size();
    const size_t nv_to_decimate =
        NV - std::max(static_cast<int>(std::lround((1 - options.strength) * NV)), 3);
    auto& vertices = mesh.vertices;
    auto& indices = mesh.indices;

    // ignoring the 4th and the following values (if exist) in vertices.
    for (auto& vert : vertices)
        vert.resize(3);

    if (nv_to_decimate == 0)
        return true;

    // keeps a reference to indices and init face2edge and edges
    Internal::Connectivity conn{indices, {}, {}};

    std::vector<bool> deleted_face(indices.size(), false);
    std::vector<bool> deleted_vertex(NV, true);

    // [1] find out information of edges (endpoints, incident faces) and face2edge
    {
        std::forward_list<idx> invalid_faces;
        Internal::construct_edges(vertices, conn, invalid_faces);

        for (auto df : invalid_faces)
            deleted_face[df] = true;
    }

    // [2] compute quadrics of vertices
    Internal::Q quadrics;
    Internal::compute_quadrics(vertices, conn, quadrics, deleted_face, options);

    // [3] assigning edge errors using quadrics
    Internal::compute_errors(vertices, quadrics, conn.edges, options.fix_boundary);

    // [4] create priority queue on quadric error
    Internal::QEMHeap heap(conn.edges, !options.fix_boundary);

    for (idx f = 0; f < indices.size(); ++f) {
        if (!deleted_face[f]) {
            const auto& face = indices[f];
            for (idx v : face)
                deleted_vertex[v] = false;
        }
    }

    size_t nv = nv_to_decimate;
    while (!heap.empty() && nv > 0) {
        // target the least-error edge, if it is what we saw last iteration,
        // it means loop should stop because all remaining edges have been penalized
        const idx target = heap.top();
        if (conn.edges[target].error >= std::numeric_limits<double>::max())
            break;

        const auto& edge = conn.edges[target];

        assert(is_valid_edge_target(edge, deleted_face, deleted_vertex,
                                    options.fix_boundary));

        // [5] collapse the least-error edge until mesh is simplified enough
        bool collapsed = edge_collapse(vertices, conn, quadrics, heap, target, options);
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
    Internal::compact_data(deleted_vertex, deleted_face, vertices, indices);

    return true;
}

bool simplify(const TriMesh& input, TriMesh& output, const SimplifyOptions& options) {
    output = input;
    return simplify(output, options);
}

} // namespace MeshSimpl
