#ifndef MESH_SIMPL_PROC_HPP
#define MESH_SIMPL_PROC_HPP

#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Edge;
class Faces;
class QEMHeap;
class Vertices;

void weightQuadrics(Quadric& quadric, double faceArea, WEIGHTING strategy);

// Compute quadrics Q for every vertex
// Output `quadrics`: one quadric for each vertex
void computeQuadrics(Vertices& vertices, Faces& faces, Edges& edges,
                     const SimplifyOptions& options);

bool edgeTopoCorrectness(const Faces& faces, const Edges& edges);

// Build connectivity from `conn.indices`
// Output `conn.edges`: list of `Edge`
// Output `conn.face2edge`: |F|x3 matrix,
//  face2edge[face index][order of vertex v] = index of edge across from v on
//  face
// Output `invalid_faces`: list of invalid face indexes (incident to
// non-manifold edges)
//
// This method does not initialize members optimal_pos and error in struct Edge.
void buildConnectivity(Vertices& vertices, Faces& faces, Edges& edges);

// Returns true if the movement of vertex will cause this face to flip too much
// to accept
bool isFaceFolded(const Vertices& vertices, const Faces& faces, idx f,
                  order moved, const vec3d& position, double angle);

bool isFaceElongated(const vec3d& pos0, const vec3d& pos1, const vec3d& pos2,
                     double ratio);

// Returns true if edge is collapsed
int edgeCollapse(Vertices& vertices, Faces& faces, QEMHeap& heap, Edge& target,
                 const SimplifyOptions& options);

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_PROC_HPP
