#ifndef LIB_MESH_SIMPL_ECOL_HPP
#define LIB_MESH_SIMPL_ECOL_HPP

#include "edge.hpp"
#include "neighbor.hpp"
#include "qem_heap.hpp"
#include "types.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

// Replace a vertex in edge, and update other members then fix priority in heap
void update_error_and_center(const V& vertices, const Q& quadrics,
                             QEMHeap& heap, Edge* edge_ptr, bool fix_boundary);

// Returns true if the movement of vertex will cause this face to flip too much
// to accept
bool face_fold_over(const V& vertices, const F& indices, const Neighbor& nb,
                    idx v_move, const vec3d& move_to, double angle);

// Returns true if the face is extremely elongated. It calculates inverse of
// triangle aspect ratio and compare to a limit
bool extremely_elongated(const V& vertices, const F& indices,
                         const Neighbor& nb, const vec3d& center_pos,
                         double ratio);

// Find relevant faces (neighbor faces of two endpoints) one by one around the
// collapsed edge. Meanwhile this function does geometry and connectivity check
// to avoid fold-over faces and non-manifold edges.
//
// Traverse neighbors around deleted vertex first, then around kept vertex,
// both clockwise. One exception is when traversing neighbors of boundary kept
// vertex. Also refer to data structure "Neighbor" in neighbor.h
//
// Outputs: v_del_neighbors -- neighbors of deleted vertex
//          c_kept_neighbor_edges -- index of edges incident to kept vertex
// Returns: true if no problem is found -- this edge collapse operation is
// acceptable
bool scan_neighbors(const V& vertices, const Connectivity& conn,
                    const Edge& edge, std::vector<Neighbor>& v_del_neighbors,
                    std::vector<idx>& v_kept_neighbor_edges, order del_ord,
                    const SimplifyOptions& options);

// Returns true if edge is collapsed
bool edge_collapse(V& vertices, Internal::Connectivity& conn, Q& quadrics,
                   QEMHeap& heap, idx ecol_target, order del_ord,
                   const SimplifyOptions& options);

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // LIB_MESH_SIMPL_ECOL_HPP
