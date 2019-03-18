//
// Created by nickl on 01/08/19.
//

#ifndef LIB_MESH_SIMPL_TYPES_H
#define LIB_MESH_SIMPL_TYPES_H

#include <array>
#include <vector>

namespace MeshSimpl {

typedef unsigned int idx;                   // edge, vertex, face index
typedef signed char order;                  // edge/vertex local order to face; [0, 3)
typedef std::array<double, 3> vec3d;        // double
typedef std::array<idx, 3> vec3i;           // idx
typedef std::array<idx, 2> vec2i;           // idx
typedef std::vector<std::vector<double>> V; // input/output vertex positions
typedef std::vector<std::vector<idx>> F;    // input/output face indices

enum WEIGHTING { UNIFORM, BY_AREA };

namespace Internal {

static const order INVALID = -1;

// A quadric Q consists of a symmetric 3x3 matrix A, a vec3 b, and a scalar c
typedef std::array<double, 10> Quadric;

// Edge defines the struct of an edge
struct Edge {
    enum BOUNDARY_V { V0 = 0, V1 = 1, BOTH, NONE };

    // sum of quadrics of two endpoints
    Quadric q;
    // where this edge collapse into
    vec3d center;

    // index of two endpoints,
    // when constructing a set of edges, this is used as key by forcing v0 < v1,
    // the increasing order is unnecessary if not in a set, but mind "boundary_v"
    vec2i vertices;

    // index of two incident faces,
    // the order of two faces is insignificant in this program (swappable)
    vec2i faces;
    // order of this edge in two incident faces,
    // match the order to member "faces"
    std::array<order, 2> idx_in_face;

    // quadric error value, insignificant if edge on boundary (boundary_v = BOTH)
    double error;

    // vertices on boundary
    BOUNDARY_V boundary_v;
};

// Used to constraint boundary from being reshaped. Necessary if boundary is not fixed
struct ConstraintPlane {
    std::vector<bool> on_boundary;
    std::vector<order> boundary_e_order;
    bool enabled;
};

typedef std::vector<Edge> E;
typedef std::vector<Quadric> Q;
typedef std::vector<vec3i> F2E;

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_TYPES_H
