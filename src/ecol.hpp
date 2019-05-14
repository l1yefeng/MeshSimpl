#ifndef MESH_SIMPL_ECOL_HPP
#define MESH_SIMPL_ECOL_HPP

#include "neighbor.hpp"
#include "qem_heap.hpp"
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

// Replace a vertex in edge, and update other members then fix priority in heap
void update_error_and_center(const V &vertices, const Q &quadrics,
                             QEMHeap &heap, Edge *edge_ptr, bool fix_boundary);

// Returns true if the movement of vertex will cause this face to flip too much
// to accept
bool face_fold_over(const V &vertices, const F &faces, const Neighbor &nb,
                    idx v_move, const vec3d &move_to, double angle);

// Returns true if the face is extremely elongated. It calculates inverse of
// triangle aspect ratio and compare to a limit
bool extremely_elongated(const V &vertices, const F &faces, const Neighbor &nb,
                         const vec3d &center_pos, double ratio);

// Returns true if edge is collapsed
bool edge_collapse(V &vertices, F &faces, Q &quadrics, QEMHeap &heap,
                   Edge &target, const SimplifyOptions &options);

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_ECOL_HPP
