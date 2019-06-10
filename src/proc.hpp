#ifndef MESH_SIMPL_PROC_HPP
#define MESH_SIMPL_PROC_HPP

#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Faces;
class Vertices;

// Weight the quadric according to strategy
void weightQuadrics(Quadric& quadric, double faceArea, WEIGHTING strategy);

// Compute quadrics Q for every vertex
void computeQuadrics(Vertices& vertices, const Faces& faces,
                     const SimplifyOptions& options);

bool edgeTopoCorrectness(const Faces& faces, const Edges& edges);

// Build connectivity, namely creating edges, assign sides to faces, set
// vertices boundary or not. This allows traversal on the mesh via:
//  * Faces::side()
//  * Faces::v()
//  * Edge::endpoint()
//  * Edge::face()
//  * Edge::ordInF()
void buildConnectivity(Vertices& vertices, Faces& faces, Edges& edges);

// Returns true if the movement of vertex will cause this face to flip too much
bool isFaceFlipped(const Vertices& vertices, const Faces& faces, idx f,
                   order moved, const vec3d& position, double angle);

// Returns true if the face formed with input positions is extremely elongated
bool isFaceElongated(const vec3d& pos0, const vec3d& pos1, const vec3d& pos2,
                     double ratio);

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_PROC_HPP
