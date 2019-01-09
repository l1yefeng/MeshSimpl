//
// Created by nickl on 01/08/19.
//

#ifndef LIB_MESH_SIMPL_TYPES_H
#define LIB_MESH_SIMPL_TYPES_H

#include <array>
#include <vector>

namespace MeshSimpl {

typedef unsigned int idx;
typedef std::array<double, 3> vec3d;        // double
typedef std::array<idx, 3> vec3i;           // idx
typedef std::array<idx, 2> vec2i;           // idx
typedef std::vector<std::vector<double>> V; // input/output vertex positions
typedef std::vector<std::vector<idx>> F;    // input/output face indices

namespace Internal {

// A quadric Q consists of a symmetric 3x3 matrix A, a vec3 b, and a scalar c
typedef std::array<double, 10> Quadric;

enum BOUNDARY_V { NONE, BOTH, V0, V1 };

// Edge defines the struct of an edge
struct Edge {
    vec2i vertices;        // index of two end vertices, unique, v0 < v1
    vec2i faces;           // index of two incident faces
    vec2i idx_in_face;     // index of this (0, 1, 2) in faces
    BOUNDARY_V boundary_v; // vertices on boundary
    vec3d center;          // where this edge collapse into
    double error;          // quadric error value, undefined if on boundary
    Quadric q;             // sum of quadrics of two vertices
};

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_TYPES_H
