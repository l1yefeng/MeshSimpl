//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_PRE_PROC_HPP
#define LIB_MESH_SIMPL_PRE_PROC_HPP

#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

void weight_quadric(Quadric& quadric, double face_area, WEIGHTING strategy);

// Compute quadrics Q for every vertex
Q compute_quadrics(const V& vertices, Internal::Connectivity& conn,
                   const SimplifyOptions& options);

// Returns face2edge and edges.
// face2edge is |F|x3 with each value indexing a unique edge in edges.
// This method does not initialize members optimal_pos and error in struct Edge.
void construct_edges(size_t vertex_cnt, Internal::Connectivity& conn);

} // namespace Internal
} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_PRE_PROC_HPP
