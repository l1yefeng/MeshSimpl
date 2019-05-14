//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_PRE_PROC_HPP
#define MESH_SIMPL_PRE_PROC_HPP

#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

void weight_quadric(Quadric &quadric, double face_area, WEIGHTING strategy);

// Compute quadrics Q for every vertex
// Output `quadrics`: one quadric for each vertex
void compute_quadrics(const V &vertices, Q &quadrics, F &faces, E &edges,
                      const SimplifyOptions &options);

// Build connectivity from `conn.indices`
// Output `conn.edges`: list of `Edge`
// Output `conn.face2edge`: |F|x3 matrix,
//  face2edge[face index][order of vertex v] = index of edge across from v on
//  face
// Output `invalid_faces`: list of invalid face indexes (incident to
// non-manifold edges)
//
// This method does not initialize members optimal_pos and error in struct Edge.
void construct_edges(const V &vertices, F &faces, E &edges);

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_PRE_PROC_HPP
