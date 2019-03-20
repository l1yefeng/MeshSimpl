//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_PRE_PROC_H
#define LIB_MESH_SIMPL_PRE_PROC_H

#include "types.hpp"

namespace MeshSimpl {

namespace Internal {

static const double CONSTRAINT_PLANE_C = 1.0;

// Compute quadrics Q for every vertex
Q compute_quadrics(const V& vertices, const F& indices,
                   const std::vector<char>& boundary_flags, WEIGHTING weighting);

// Returns face2edge and edges.
// face2edge is |F|x3 with each value indexing a unique edge in edges.
// This method does not initialize members optimal_pos and error in struct Edge.
std::pair<E, std::vector<vec3i>> construct_edges(const F& indices, size_t vertex_cnt,
                                                 std::vector<char>& boundary_flags);

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_PRE_PROC_H
