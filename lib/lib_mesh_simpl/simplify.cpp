//
// Created by nickl on 1/8/19.
//

#include "simplify.h"
#include "ecol.h"
#include "measure.h"
#include "post_proc.h"
#include "pre_proc.h"
#include <iostream>
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

    Internal::Measure measure;
    measure.start("all");
    if (options.debug)
        std::cout << "INFO::INPUT: mesh contains " << vertices.size() << " vertices, "
                  << indices.size() << " faces" << std::endl;

    const size_t NV = vertices.size();
    const size_t nv_to_decimate =
        NV - std::max(static_cast<size_t>(std::lround((1 - options.strength) * NV)),
                      Internal::MIN_NR_VERTICES);
    auto out_vertices = vertices;
    auto out_indices = indices;
    if (options.debug)
        std::cout << "INFO::INPUT: simplification strength is " << options.strength
                  << ", process will stop when the number of remaining vertices is "
                  << (NV - nv_to_decimate) << ", or when no remaining edge can be contracted"
                  << std::endl;

    if (nv_to_decimate == 0)
        return {out_vertices, out_indices};

    measure.start("compute_quadrics");
    auto quadrics = Internal::compute_quadrics(vertices, indices, options.weight_by_face);
    auto ms = measure.stop("compute_quadrics");
    if (options.debug)
        std::cout << "INFO::PROCESSING::PRE_PROCESS: computing vertex quadrics completed (" << ms
                  << " ms)" << std::endl;

    measure.start("construct_edges");
    auto edge_topo = Internal::construct_edges(indices, NV);
    ms = measure.stop("construct_edges");
    if (options.debug)
        std::cout << "INFO::PROCESSING::PRE_PROCESS: computing connectivity information completed ("
                  << ms << " ms)" << std::endl;
    std::vector<Internal::Edge>& edges = edge_topo.first;
    std::vector<vec3i>& face2edge = edge_topo.second;

    measure.start("heap");
    Internal::compute_errors(vertices, quadrics, edges);
    // create priority queue on quadric error
    Internal::QEMHeap heap(edges);
    ms = measure.stop("heap");
    if (options.debug)
        std::cout << "INFO::PROCESSING: building binomial min-heap of edges keying on quadric "
                     "error completed ("
                  << ms << " ms)" << std::endl;

    std::vector<bool> deleted_vertex(NV, false);
    std::vector<bool> deleted_face(indices.size(), false);

    // one run of a series of edge collapse that is supposed to decimate nv_decimate vertices and
    // returns true if it ends because this target is achieved instead of because of other reason
    size_t counter = 0;
    size_t nv = nv_to_decimate;
    if (options.debug)
        std::cout << "INFO::PROCESSING: starting edge-collapse iterations ..." << std::endl;
    measure.start("edge_collapse");
    while (!heap.empty() && nv > 0) {
        if (options.debug && (++counter) % (nv_to_decimate / 5) == 0)
            std::cout << "INFO::PROCESSING::PROGRESS: running edge-collapse operator -- "
                      << (counter / (nv_to_decimate / 5) * 20) << "%" << std::endl;

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
            --nv;
        }
    }
    ms = measure.stop("edge_collapse");
    if (options.debug)
        std::cout << "INFO::PROCESSING: edge-collapse iterations completed (" << ms << " ms)"
                  << std::endl;

    // edges are pointless from this point on, but need to fix vertices and indices
    measure.start("compact_data");
    Internal::compact_data(deleted_vertex, deleted_face, out_vertices, out_indices);
    ms = measure.stop("compact_data");
    if (options.debug)
        std::cout
            << "INFO::PROCESSING::POST_PROCESS: compacting vertices and indices data completed ("
            << ms << " ms)" << std::endl;

    ms = measure.stop("all");
    if (options.debug)
        std::cout << "INFO::FINISH: simplification ended -- output mesh contains "
                  << out_vertices.size() << " vertices, " << out_indices.size() << " faces (" << ms
                  << " ms)" << std::endl;
    return {out_vertices, out_indices};
}

} // namespace MeshSimpl
