#ifndef LIB_MESH_SIMPL_ECOL_HPP
#define LIB_MESH_SIMPL_ECOL_HPP

#include "edge.hpp"
#include "neighbor.hpp"
#include "qem_heap.hpp"
#include "types.hpp"
#include "util.hpp"

namespace MeshSimpl {

namespace Internal {

static const double FOLD_OVER_COS_ANGLE = std::cos(160);
static const size_t ESTIMATE_VALENCE = 8;
static const double ASPECT_RATIO_AT_LEAST = 0.01;

// Set error and center of edge by choosing a position to collapse into
void optimal_ecol_vertex_placement(const V& vertices, Edge& edge);

// Set quadric, error, and collapse center for a given edge
void set_edge_error(const V& vertices, const Q& quadrics, Edge& edge, bool fix_boundary);

// Compute quadric error for every edge at initialization
void compute_errors(const V& vertices, const Q& quadrics, E& edges, bool fix_boundary);

// Replace a vertex in edge, and update other members then fix priority in heap
void update_error_and_center(const V& vertices, const Q& quadrics, QEMHeap& heap,
                             Edge* const edge_ptr, bool fix_boundary);

// Returns true if the movement of vertex will cause this face to flip too much to accept
bool face_fold_over(const V& vertices, const F& indices, const Neighbor& nb, idx v_move,
                    const vec3d& move_to);

// Returns true if the face is extremely elongated. It calculates inverse of
// triangle aspect ratio and compare to a limit
bool extremely_elongated(const V& vertices, const F& indices, const Neighbor& nb,
                         const vec3d& center_pos);

bool boundary_scan_neighbors(const V& vertices, const F& indices, const E& edges,
                             const F2E& face2edge, const Edge& edge,
                             std::vector<Neighbor>& v_del_neighbors,
                             std::vector<idx>& v_kept_neighbor_edges);

bool boundary_edge_collapse(V& vertices, F& indices, E& edges, F2E& face2edge,
                            Q& quadrics, QEMHeap& heap, const idx ecol_target);

// Find relevant faces (neighbor faces of two endpoints) one by one around the
// collapsed edge. Meanwhile this function does geometry and connectivity check to
// avoid fold-over faces and non-manifold edges.
//
// Traverse neighbors around deleted vertex first, then around kept vertex,
// both clockwise. One exception is when traversing neighbors of boundary kept vertex.
// Also refer to data structure "Neighbor" in neighbor.h
//
// Outputs: v_del_neighbors -- neighbors of deleted vertex
//          c_kept_neighbor_edges -- index of edges incident to kept vertex
// Returns: true if no problem is found -- this edge collapse operation is acceptable
bool scan_neighbors(const V& vertices, const F& indices, const E& edges,
                    const F2E& face2edge, const Edge& edge,
                    std::vector<Neighbor>& v_del_neighbors,
                    std::vector<idx>& v_kept_neighbor_edges);

// Returns true if edge is collapsed
bool edge_collapse(V& vertices, F& indices, E& edges, F2E& face2edge, Q& quadrics,
                   QEMHeap& heap, const idx ecol_target, bool fix_boundary);

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_ECOL_HPP
