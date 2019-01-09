//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_PRE_PROC_H
#define LIB_MESH_SIMPL_PRE_PROC_H

#include "types.h"

namespace MeshSimpl
{
namespace Internal
{

void compute_quadrics_per_face(const V& vertices, const std::vector<idx>& face,
                               bool weight_by_area, std::vector<Quadric>& quadrics);

// Compute quadrics Q for every vertex. If not weighted by area then it's uniform
std::vector<Quadric>
compute_quadrics(const V& vertices, const F& indices, bool weight_by_area);

// Compute quadrics Q for every vertex. If not weighted by area then it's uniform
std::vector<Quadric>
compute_quadrics(const V& vertices, const F& indices, const std::vector<bool>& deleted_face,
                 bool weight_by_area);

// Returns face2edge and edges. face2edge is NFx3 with each value indexing a unique edge in edges.
// This method does not initialize members optimal_pos and error in struct Edge.
std::pair<std::vector<Edge>, std::vector<vec3i>>
construct_edges(const F& indices, size_t vertex_cnt);

}
}

#endif // LIB_MESH_SIMPL_PRE_PROC_H
