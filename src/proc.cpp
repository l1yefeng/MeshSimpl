//
// Created by nickl on 1/8/19.
//

#include <array>
#include <cassert>
#include <initializer_list>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "quadric.hpp"
#include "edge.hpp"
#include "faces.hpp"
#include "proc.hpp"
#include "util.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

void computeQuadrics(Vertices &vertices, const Faces &faces,
                     const SimplifyOptions &options) {
  for (idx f = 0; f < faces.size(); ++f) {
    // calculate the plane of this face (n and d: n'v+d=0 defines the plane)
    const vec3d edgeVec2 = faces.edgeVec(f, 2, vertices);
    const vec3d edgeVec1 = faces.edgeVec(f, 1, vertices);
    vec3d normal = cross(edgeVec2, edgeVec1);
    // |normal| = area, used for normalization and weighting quadrics
    const double area = magnitude(normal);
    if (area != 0)
      normal /= area;
    else
      continue;

    // d = -n*v0
    const double d = -dot(normal, faces.vPos(f, 0, vertices));

    // calculate quadric Q = (A, b, c) = (nn', dn, d*d)
    Quadric q (normal, d);
    if (options.weightByArea) q *= area;

    for (order k : {0, 1, 2}) vertices.increaseQ(faces.v(f, k), q);
  }

  if (!options.fixBoundary) {
    for (idx f = 0; f < faces.size(); ++f) {
      if (!faces.onBoundary(f)) continue;

      const std::array<vec3d, 3> e{faces.edgeVec(f, 0, vertices),
                                   faces.edgeVec(f, 1, vertices),
                                   faces.edgeVec(f, 2, vertices)};
      const vec3d nFace = cross(e[0], e[1]);

      for (order k : {0, 1, 2}) {
        if (!faces.side(f, k)->onBoundary()) continue;

        vec3d normal = cross(nFace, e[k]);
        const double normalMag = magnitude(normal);
        if (normalMag != 0)
          normal /= normalMag;
        else
          continue;

        const double d = -dot(normal, faces.vPos(f, next(k), vertices));
        Quadric q(normal, d);
        q *= options.borderConstraint;
        if (options.weightByArea) q *= magnitude(nFace);

        vertices.increaseQ(faces.v(f, next(k)), q);
        vertices.increaseQ(faces.v(f, prev(k)), q);
      }
    }
  }
}

bool edgeTopoCorrectness(const Faces &faces, const Edges &edges) {
  for (idx f = 0; f < faces.size(); ++f) {
    for (order ord = 0; ord < 3; ++ord) {
      auto &vv = faces.side(f, ord)->endpoints();
      assert(vv[0] < vv[1]);
      idx vSmaller = faces.v(f, next(ord));
      idx vLarger = faces.v(f, prev(ord));
      if (vSmaller > vLarger) std::swap(vSmaller, vLarger);

      // edge vertices should match face corners
      assert(vv[0] == vSmaller && vv[1] == vLarger);
    }
  }

  return true;
}

void buildConnectivity(Vertices &vertices, Faces &faces, Edges &edges) {
  std::map<std::pair<idx, idx>, std::pair<Edge, bool>> edgeSet;

  // insert all edges into edgeSet and find out if it is on boundary
  for (idx f = 0; f < faces.size(); ++f) {
    const auto &face = faces[f];
    for (order k = 0; k < 3; ++k) {
      // construct edge (v[i], v[j]);
      // edge local index will be k (= that of the 3rd vertex)
      const order i = next(k);
      const order j = prev(k);
      idx v0 = face[i];
      idx v1 = face[j];
      if (v0 > v1) std::swap(v0, v1);
      auto res = edgeSet.emplace(std::make_pair(v0, v1),
                                 std::make_pair(Edge(vertices, v0, v1), false));
      auto inserted = res.second;
      auto &edge = res.first->second.first;
      auto &twoWings = res.first->second.second;

      if (inserted) {
        edge.setWing(0, f, k);
        assert(!twoWings);
      } else {
        if (twoWings) {
          std::stringstream ss;
          ss << "ERROR::INPUT_MESH: found non-manifold edge" << std::endl;
          for (idx _f : {edge.face(0), edge.face(1), f}) {
            ss << "                   face #" << _f << ":";
            for (int _i = 0; _i < 3; ++_i) ss << "\t" << faces[_f][_i] + 1;
            ss << std::endl;
          }
          throw std::invalid_argument(ss.str());
        } else {
          edge.setWing(1, f, k);
          twoWings = true;
        }
      }
    }
  }

  // populate edges vector from map
  edges.reserve(edgeSet.size());
  for (const auto &elem : edgeSet) {
    auto &edge = elem.second.first;
    auto &twoWings = elem.second.second;

    edges.emplace_back(edge);

    faces.setSide(edge.face(0), edge.ordInF(0), &edges.back());
    if (!twoWings) {
      for (order i : {0, 1}) vertices.setBoundary(edge.endpoint(i), true);
    } else {
      faces.setSide(edge.face(1), edge.ordInF(1), &edges.back());
    }
  }

  assert(edgeTopoCorrectness(faces, edges));
}

bool isFaceFlipped(const Vertices &vertices, const Faces &faces, idx f,
                   order moved, const vec3d &position, double angle) {
  const vec3d &vk = faces.vPos(f, moved, vertices);
  const vec3d &vi = faces.vPos(f, next(moved), vertices);
  const vec3d &vj = faces.vPos(f, prev(moved), vertices);
  const vec3d edgeVec0 = vj - vi;
  const vec3d edgeVec1 = vk - vj;
  const vec3d edgeVec1New = position - vj;
  vec3d normalPrv = cross(edgeVec0, edgeVec1);
  double magPrv = magnitude(normalPrv);
  if (magPrv != 0) {
    vec3d normalNew = cross(edgeVec0, edgeVec1New);
    double magNew = magnitude(normalNew);
    if (magNew == 0) return true;
    normalPrv /= magPrv;
    normalNew /= magNew;
    double cos = dot(normalPrv, normalNew);
    return cos < angle;
  } else {
    return false;
  }
}

}  // namespace Internal
}  // namespace MeshSimpl
