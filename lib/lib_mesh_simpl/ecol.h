#ifndef LIB_MESH_SIMPL_ECOL_H
#define LIB_MESH_SIMPL_ECOL_H

#include "qem_heap.h"
#include "types.h"
#include "util.h"
#include <queue>

namespace MeshSimpl {

namespace Internal {

static const size_t MIN_NR_VERTICES = 4;
static const double FOLD_OVER_COS_ANGLE = std::cos(160);

// Set error and center of edge by choosing a position to collapse into
void optimal_ecol_vertex_placement(const V& vertices, Edge& edge);

// Set quadric, error, and collapse center for a given edge
void set_edge_error(const V& vertices, const std::vector<Quadric>& quadrics, Edge& edge);

// Compute quadric error for every edge at initialization
void compute_errors(const V& vertices, const std::vector<Quadric>& quadrics,
                    std::vector<Edge>& edges);

// Recompute quadric error for evert edge with idx in range
void recompute_errors(const V& vertices, const std::vector<Quadric>& quadrics,
                      std::vector<Edge>& edges, std::vector<idx>::const_iterator edge_idx_begin,
                      std::vector<idx>::const_iterator edge_idx_end);

// Returns true if the movement of vertex will cause this face to flip too much to accept
bool face_fold_over(const V& vertices, idx v0, idx v1, idx v2_prev, const vec3d& v2_new_pos);

// Replace a vertex in edge, and update other members then fix priority in heap
void update_error_and_center(const V& vertices, const std::vector<Quadric>& quadrics, QEMHeap& heap,
                             Edge* edge_ptr);

// Find relevant faces one by one around the collapsed edge and push them into two queues.
// For each face around v_del/v_kept, a 3-tuple called `fve' is stored in the corresponding queue;
//      `f' - idx of face;
//      `e' - local idx (0/1/2) of the edge to this face, shared by this face and previous face';
//      `v' - local idx (0/1/2) of the vertex to this face, on edge `e' and not centered vertex;
bool scan_neighbors(const V& vertices, const F& indices, const std::vector<Edge>& edges,
                    const std::vector<vec3i>& face2edge, const Edge& edge, idx v_del, idx v_kept,
                    std::queue<vec3i>& fve_queue_v_del, std::queue<vec3i>& fve_queue_v_kept);

// Returns true if topology will be preserved
bool topology_preserved(const std::vector<vec3i>& face2edge,
                        const std::queue<vec3i>& fve_queue_v_del,
                        const std::queue<vec3i>& fve_queue_v_kept);

inline idx vi_in_face(const F& indices, const idx f, const idx v) {
    if (indices[f][0] == v)
        return 0;
    if (indices[f][1] == v)
        return 1;
    assert(indices[f][2] == v);
    return 2;
}

inline idx fi_in_edge(const Edge& edge, const idx f) {
    if (edge.faces[0] == f)
        return 0;
    assert(edge.faces[1] == f);
    return 1;
};

inline idx fve_center(const vec3i& fve) { return 3 - fve[1] - fve[2]; }

// Returns true if edge is collapsed
bool collapse_interior_edge(V& vertices, F& indices, std::vector<Edge>& edges,
                            std::vector<vec3i>& face2edge, std::vector<Quadric>& quadrics,
                            std::vector<bool>& deleted_vertex, std::vector<bool>& deleted_face,
                            QEMHeap& heap, idx ecol_target);

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_ECOL_H
