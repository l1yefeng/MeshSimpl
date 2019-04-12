//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_PRE_PROC_HPP
#define LIB_MESH_SIMPL_PRE_PROC_HPP

#include "types.hpp"
#include <forward_list>
#include <map>

namespace MeshSimpl {
namespace Internal {

void weight_quadric(Quadric& quadric, double face_area, WEIGHTING strategy);

// Compute quadrics Q for every vertex
// Output `quadrics`: one quadric for each vertex
void compute_quadrics(const V& vertices, Internal::Connectivity& conn, Q& quadrics,
                      const std::vector<bool>& deleted_face,
                      const SimplifyOptions& options);

// Short-hand type definition
typedef std::forward_list<std::pair<idx, order>> flist_t;

// Ensures there is no edge with more than two incident faces.
// Output `invalid_faces`: list of invalid face indexes (incident to non-manifold edges)
//
// Afterwards, it is possible some edge has two INVALID `ord_in_face` value
// because its incident faces were made invalid.
//
// The definition looks awfully complex but it is simply removing large-area faces
// incident on non-manifold edges until two are left. Most codes are updating
// connectivity data structures.
void ensure_manifold(std::map<std::pair<idx, idx>, std::pair<Edge, flist_t>>& edge_set,
                     const V& vertices, const F& indices,
                     std::forward_list<idx>& invalid_faces);

// Build connectivity from `conn.indices`
// Output `conn.edges`: list of `Edge`
// Output `conn.face2edge`: |F|x3 matrix,
//  face2edge[face index][order of vertex v] = index of edge across from v on face
// Output `invalid_faces`: list of invalid face indexes (incident to non-manifold edges)
//
// This method does not initialize members optimal_pos and error in struct Edge.
void construct_edges(const V& vertices, Internal::Connectivity& conn,
                     std::forward_list<idx>& invalid_faces);

} // namespace Internal
} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_PRE_PROC_HPP
