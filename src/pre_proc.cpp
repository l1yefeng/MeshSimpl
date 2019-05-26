//
// Created by nickl on 1/8/19.
//

#include "pre_proc.hpp"
#include <array>             // for array
#include <cassert>           // for assert
#include <initializer_list>  // for initializer_list
#include <limits>            // for numeric_limits
#include <map>               // for map, _Rb_tree_iterator
#include <sstream>           // for operator<<, basic_ostream, stringstream
#include <stdexcept>         // for invalid_argument
#include <utility>           // for pair, swap, make_pair
#include <vector>            // for vector, vector<>::reference, _Bit_reference
#include "edge.hpp"          // for Edge
#include "faces.hpp"         // for Faces
#include "util.hpp"          // for operator*=, next, cross, magnitude, prev
#include "vertices.hpp"      // for Vertices

namespace MeshSimpl {
namespace Internal {

void weight_quadric(Quadric &quadric, double face_area, WEIGHTING strategy) {
  assert(face_area != 0);
  switch (strategy) {
    case BY_AREA:
      quadric *= face_area;
      break;
    case BY_AREA_INV:
      if (face_area <= std::numeric_limits<double>::epsilon())
        quadric *= 0.0;
      else
        quadric *= 1 / face_area;
      break;
    case UNIFORM:
    default:
      if (face_area <= std::numeric_limits<double>::epsilon()) quadric *= 0.0;
      break;
  }
}

void compute_quadrics(Vertices &vertices, Faces &faces, E &edges,
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
    Quadric q = make_quadric(normal, d);
    weight_quadric(q, area, options.weighting);

    for (order k : {0, 1, 2}) vertices.increaseQ(faces.v(f, k), q);
  }

  if (!options.fix_boundary) {
    for (idx f = 0; f < faces.size(); ++f) {
      if (!faces.onBoundary(f)) continue;

      const std::array<vec3d, 3> e{faces.edgeVec(f, 0, vertices),
                                   faces.edgeVec(f, 1, vertices),
                                   faces.edgeVec(f, 2, vertices)};
      const vec3d n_face = cross(e[0], e[1]);

      for (order k : {0, 1, 2}) {
        if (!faces.side(f, k)->on_boundary()) continue;

        vec3d normal = cross(n_face, e[k]);
        const double normal_mag = magnitude(normal);
        if (normal_mag != 0)
          normal /= normal_mag;
        else
          continue;

        const double d = -dot(normal, faces.vPos(f, next(k), vertices));
        Quadric q = make_quadric(normal, d);
        q *= options.border_constraint;
        weight_quadric(q, magnitude(n_face), options.weighting);

        vertices.increaseQ(faces.v(f, next(k)), q);
        vertices.increaseQ(faces.v(f, prev(k)), q);
      }
    }
  }
}

bool edge_topo_correctness(const Faces &faces, const E &edges) {
  for (idx f = 0; f < faces.size(); ++f) {
    for (order ord = 0; ord < 3; ++ord) {
      auto &vv = faces.side(f, ord)->vertices();
      assert(vv[0] < vv[1]);
      idx v_smaller = faces.v(f, next(ord));
      idx v_larger = faces.v(f, prev(ord));
      if (v_smaller > v_larger) std::swap(v_smaller, v_larger);

      // edge vertices should match face corners
      assert(vv[0] == v_smaller && vv[1] == v_larger);
    }
  }

  return true;
}

void construct_edges(Vertices &vertices, Faces &faces, E &edges) {
  std::map<std::pair<idx, idx>, Edge> edge_set;

  // insert all edges into edge_set and find out if it is on boundary
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
      auto res = edge_set.emplace(std::make_pair(v0, v1), Edge(v0, v1));
      auto inserted = res.second;
      auto &edge = res.first->second;

      if (inserted) {
        edge.attach_1st_face(f, k);
      } else {
        bool ok = edge.attach_2nd_face(f, k);
        if (!ok) {
          std::stringstream ss;
          ss << "ERROR::INPUT_MESH: found non-manifold edge" << std::endl;
          for (idx _f : {edge.face(0), edge.face(1), f}) {
            ss << "                   face #" << _f << ":";
            for (int _i = 0; _i < 3; ++_i) ss << "\t" << faces[_f][_i] + 1;
            ss << std::endl;
          }
          throw std::invalid_argument(ss.str());
        }
      }
    }
  }

  // populate edges vector from map
  edges.reserve(edge_set.size());
  for (const auto &elem : edge_set) edges.emplace_back(elem.second);

  std::vector<bool> vertex_on_boundary(vertices.size(), false);

  for (auto &edge : edges) {
    // identify boundary vertices
    if (edge.both_v_on_border()) {
      vertex_on_boundary[edge[0]] = true;
      vertex_on_boundary[edge[1]] = true;
    }

    // populate face2edge references
    faces.setSide(edge.face(0), edge.ord_in_face(0), &edge);
    if (!edge.both_v_on_border())
      faces.setSide(edge.face(1), edge.ord_in_face(1), &edge);
  }

  // non-boundary edges may have one vertex on boundary, find them in this loop
  for (auto &edge : edges) {
    edge.set_v_on_border(vertex_on_boundary[edge[0]],
                         vertex_on_boundary[edge[1]]);
  }

  assert(edge_topo_correctness(faces, edges));
}

}  // namespace Internal
}  // namespace MeshSimpl
