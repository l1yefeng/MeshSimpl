#ifndef MESH_SIMPL_ECOL_HPP
#define MESH_SIMPL_ECOL_HPP

#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Edge;
class Faces;
class QEMHeap;
class Vertices;

// Replace a vertex in edge, and update other members then fix priority in heap
void update_error_and_center(Vertices& vertices, QEMHeap& heap, Edge* edge_ptr,
                             bool fix_boundary);

// Returns true if the movement of vertex will cause this face to flip too much
// to accept
bool is_face_folded(const Vertices& vertices, const Faces& faces, idx f,
                    order moved, const vec3d& position, double angle);

bool is_face_elongated(const vec3d& pos0, const vec3d& pos1, const vec3d& pos2,
                       double ratio);

// Returns true if edge is collapsed
int edge_collapse(Vertices& vertices, Faces& faces, QEMHeap& heap, Edge& target,
                  const SimplifyOptions& options);

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_ECOL_HPP
