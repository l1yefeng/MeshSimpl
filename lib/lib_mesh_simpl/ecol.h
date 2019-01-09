#ifndef LIB_MESH_SIMPL_ECOL_H
#define LIB_MESH_SIMPL_ECOL_H

#include "types.h"
#include "util.h"
#include "qem_heap.h"
#include <set>
#include <queue>
#include <limits>
#include <functional>

namespace MeshSimpl
{
namespace Internal
{

static const size_t MIN_NR_VERTICES = 4;
static const double FOLD_OVER_COS_ANGLE = std::cos(160);

// Set error and center of edge by choosing a position to collapse into
void optimal_ecol_vertex_placement(const V& vertices, Edge& edge);

void set_edge_error(const V& vertices, const std::vector<Quadric>& quadrics, Edge& edge);

// Compute quadric error for every edge at initialization
void init_edge_errors(const V& vertices,
                      const std::vector<Quadric>& quadrics,
                      std::vector<Edge>& edges);

bool face_fold_over(const V& vertices, idx v0, idx v1, idx v2_prev,
                    const vec3d& v2_new_pos);

// Replace a vertex in edge, and update other members then fix priority in heap
void update_error_and_center(const V& vertices, const std::vector<Quadric>& quadrics,
                             QEMHeap& heap, Edge* edge_ptr);

bool topology_preservation_check(const std::vector<vec3i>& face2edge,
                                 const std::queue<vec3i>& fve_queue_v_del,
                                 const std::queue<vec3i>& fve_queue_v_kept);

inline idx vi_in_face(const F& indices, const idx f, const idx v)
{
    if (indices[f][0] == v)
        return 0;
    if (indices[f][1] == v)
        return 1;
    assert(indices[f][2] == v);
    return 2;
}

inline idx fi_in_edge(const Edge& edge, const idx f)
{
    if (edge.faces[0] == f)
        return 0;
    assert(edge.faces[1] == f);
    return 1;
};

// Returns true if edge is collapsed
bool collapse_interior_edge(V& vertices, F& indices,
                            std::vector<Edge>& edges,
                            std::vector<vec3i>& face2edge,
                            std::vector<Quadric>& quadrics,
                            std::vector<bool>& deleted_vertex,
                            std::vector<bool>& deleted_face,
                            const std::function<void(vec3i&)>& iter_next,
                            QEMHeap& heap, idx ecol_target);

} // namespace MeshSimpl::Internal
} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_ECOL_H
