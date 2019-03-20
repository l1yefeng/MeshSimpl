//
// Created by nickl on 01/08/19.
//

#ifndef LIB_MESH_SIMPL_TYPES_H
#define LIB_MESH_SIMPL_TYPES_H

#include <array>
#include <vector>

namespace MeshSimpl {

typedef unsigned int idx;                   // edge, vertex, face index
typedef char order;                         // edge/vertex local order to face; [0, 3)
typedef std::array<double, 3> vec3d;        // double
typedef std::array<idx, 3> vec3i;           // idx
typedef std::array<idx, 2> vec2i;           // idx
typedef std::vector<std::vector<double>> V; // input/output vertex positions
typedef std::vector<std::vector<idx>> F;    // input/output face indices

enum WEIGHTING { UNIFORM, BY_AREA };

namespace Internal {

// A quadric Q consists of a symmetric 3x3 matrix A, a vec3 b, and a scalar c
typedef std::array<double, 10> Quadric;

// Defined in edge.hpp
struct Edge;

typedef std::vector<Edge> E;
typedef std::vector<Quadric> Q;
typedef std::vector<vec3i> F2E;

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_TYPES_H
