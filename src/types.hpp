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

struct SimplifyOptions {
  // simplifies until face count is 1-strength of the original,
  // only accept value in range [0, 1)
  float strength = 0.0f;

  // weight the quadrics by triangle area, i.e., Q becomes Q * weight
  // a larger weight makes the computed error larger thus "later" to modify
  // during the iterations of edge collapse operations.
  bool weightByArea = false;

  // when "fixBoundary" is true, we completely do not collapse anything on
  // boundary; otherwise, we add a "constraint plane" that is perpendicular to
  // boundary face to increase the quadric/error of boundary vertices ref:
  // Simplifying Surfaces with Color and Texture using Quadric Error Metrics
  bool fixBoundary = false;

  bool topologyModifiable = false;

  // the following are very fine grained configuration options

  // the constant that decides the weight of constraint planes (if fixBoundary
  // unset); larger borderConstraint -> harder for boundary to deform
  float borderConstraint = 1.0f;

  // used to check if, during edge collapse, faces get folded;
  double foldOverAngleThreshold = std::cos(80.0 / 180.0 * std::acos(-1));

  // used to check if a face is extremely elongated;
  // AR = 8*(s-a)*(s-b)*(s-c) / (a*b*c) where s = (a+b+c) / 2
  // BUT please provide a reciprocal, i.e., op will fail when
  //    AR > 1 / aspectRatioThreshold
  // if aspectRatioThreshold <= 0, this checking is skipped
  double aspectRatioThreshold = 0.02;
};

namespace Internal {

// Defined in edge.hpp
class Edge;

typedef std::vector<Edge> Edges;

}  // namespace Internal

}  // namespace MeshSimpl

#endif  // MESH_SIMPL_TYPES_HPP
