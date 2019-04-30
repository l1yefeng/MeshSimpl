//
// Created by nickl on 01/08/19.
//

#ifndef LIB_MESH_SIMPL_TYPES_HPP
#define LIB_MESH_SIMPL_TYPES_HPP

#include <array>
#include <cmath>
#include <vector>

namespace MeshSimpl {

typedef unsigned int idx;                   // edge, vertex, face index
typedef char order;                         // edge/vertex local order to face; [0, 3)
typedef std::array<double, 3> vec3d;        // double
typedef std::array<idx, 3> vec3i;           // idx
typedef std::array<idx, 2> vec2i;           // idx
typedef std::vector<std::vector<double>> V; // input/output vertex positions
typedef std::vector<std::vector<idx>> F;    // input/output face indices

enum WEIGHTING { UNIFORM, BY_AREA, BY_AREA_INV };

struct TriMesh {
    V vertices;
    F indices;
};

struct SimplifyOptions {
    // simplifies until vertex count is 1-strength of the original,
    // only accept value in range [0, 1)
    float strength = 0.5f;

    // weight the quadrics by triangle area, i.e., Q becomes Q * weight
    // a larger weight makes the computed error larger thus "later" to modify
    // during the iterations of edge collapse operations.
    //  - UNIFORM: no weighting
    //  - BY_AREA: larger face -> larger error
    //  - BY_AREA_INV: larger face -> smaller error
    WEIGHTING weighting = UNIFORM;

    // when "fix_boundary" is true, we completely do not collapse anything on boundary;
    // otherwise, we add a "constraint plane" that is perpendicular to boundary face
    // to increase the quadric/error of boundary vertices
    // ref: Simplifying Surfaces with Color and Texture using Quadric Error Metrics
    bool fix_boundary = false;

    // the following are very fine grained configuration options

    // the constant that decides the weight of constraint planes (if fix_boundary unset);
    // larger border_constraint -> harder for boundary to deform
    float border_constraint = 2.0f;

    // used to check if, during edge collapse, faces get folded;
    double fold_over_angle_threshold = std::cos(160);

    // used to check if, during edge collapse, faces become extremely elongated;
    // aspect_ratio = 8(s-a)(s-b)(s-c)/abc, faces with lower aspect ratio -> lower quality
    double aspect_ratio_at_least = 0.02;
};

namespace Internal {

// A quadric Q consists of a symmetric 3x3 matrix A, a vec3 b, and a scalar c
typedef std::array<double, 10> Quadric;

// Defined in edge.hpp
struct Edge;

// Defined in connectivity.hpp
struct Connectivity;

typedef std::vector<Edge> E;
typedef std::vector<Quadric> Q;
typedef std::vector<vec3i> F2E;

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_TYPES_HPP
