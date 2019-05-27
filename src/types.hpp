//
// Created by nickl on 01/08/19.
//

#ifndef MESH_SIMPL_TYPES_HPP
#define MESH_SIMPL_TYPES_HPP

#include <array>
#include <cmath>
#include <vector>

namespace MeshSimpl {

typedef unsigned int idx;             // edge, vertex, face index
typedef char order;                   // edge/vertex local order to face; [0, 3)
typedef std::array<double, 3> vec3d;  // double
typedef std::array<idx, 3> vec3i;     // idx
typedef std::array<idx, 2> vec2i;     // idx

typedef std::vector<vec3i> Indices;
typedef std::vector<vec3d> Positions;

static const order INVALID = -1;

enum WEIGHTING { UNIFORM, BY_AREA, BY_AREA_INV };

struct SimplifyOptions {
  // simplifies until vertex count is 1-strength of the original,
  // only accept value in range [0, 1)
  float strength = 0.0f;

  // weight the quadrics by triangle area, i.e., Q becomes Q * weight
  // a larger weight makes the computed error larger thus "later" to modify
  // during the iterations of edge collapse operations.
  //  - UNIFORM: no weighting
  //  - BY_AREA: larger face -> larger error
  //  - BY_AREA_INV: larger face -> smaller error
  WEIGHTING weighting = UNIFORM;

  // when "fixBoundary" is true, we completely do not collapse anything on
  // boundary; otherwise, we add a "constraint plane" that is perpendicular to
  // boundary face to increase the quadric/error of boundary vertices ref:
  // Simplifying Surfaces with Color and Texture using Quadric Error Metrics
  bool fixBoundary = false;

  // the following are very fine grained configuration options

  // the constant that decides the weight of constraint planes (if fixBoundary
  // unset); larger borderConstraint -> harder for boundary to deform
  float borderConstraint = 2.0f;

  // used to check if, during edge collapse, faces get folded;
  double foldOverAngleThreshold = std::cos(160);

  // used to check if, during edge collapse, faces become extremely elongated;
  // aspect ratio = 8(s-a)(s-b)(s-c)/abc, faces with lower aspect ratio -> lower
  // quality
  double aspectRatioAtLeast = 0.02;
};

namespace Internal {

// A quadric Q consists of a symmetric 3x3 matrix A, a vec3 b, and a scalar c
typedef std::array<double, 10> Quadric;

// Defined in edge.hpp
class Edge;

typedef std::vector<Edge> Edges;

}  // namespace Internal

}  // namespace MeshSimpl

#endif  // MESH_SIMPL_TYPES_HPP
