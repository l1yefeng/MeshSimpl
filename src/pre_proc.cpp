//
// Created by nickl on 1/8/19.
//

#include "pre_proc.hpp"
#include <map>
#include <sstream>
#include "face.hpp"
#include "util.hpp"

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

void compute_quadrics(const V &vertices, Q &quadrics, F &faces, E &edges,
                      const SimplifyOptions &options) {
  // quadrics are initialized with all zeros
  quadrics.resize(vertices.size());

  for (const auto &face : faces) {
    // calculate the plane of this face (n and d: n'v+d=0 defines the plane)
    const vec3d edge01 = vertices[face[1]] - vertices[face[0]];
    const vec3d edge02 = vertices[face[2]] - vertices[face[0]];
    vec3d normal = cross(edge01, edge02);
    // |normal| = area, used for normalization and weighting quadrics
    const double area = magnitude(normal);
    if (area != 0)
      normal /= area;
    else
      continue;

    // d = -n*v0
    const double d = -dot(normal, vertices[face[0]]);

    // calculate quadric Q = (A, b, c) = (nn', dn, d*d)
    Quadric q = make_quadric(normal, d);
    weight_quadric(q, area, options.weighting);

    for (const auto &v : face.vertices()) quadrics[v] += q;
  }

  if (!options.fix_boundary) {
    for (auto &face : faces) {
      order k;
      for (k = 0; k < 3; ++k)
        if (face.edge(k)->on_boundary()) break;

      if (k == 3) continue;

      const std::array<vec3d, 3> e{vertices[face[2]] - vertices[face[1]],
                                   vertices[face[0]] - vertices[face[2]],
                                   vertices[face[1]] - vertices[face[0]]};
      const vec3d n_face = cross(e[0], e[1]);

      for (k = 0; k < 3; ++k) {
        if (!face.edge(k)->on_boundary()) continue;

        vec3d normal = cross(n_face, e[k]);
        const double normal_mag = magnitude(normal);
        if (normal_mag != 0)
          normal /= normal_mag;
        else
          continue;

        const double d = -dot(normal, vertices[face[(k + 1) % 3]]);
        Quadric q = make_quadric(normal, d);
        q *= options.border_constraint;
        weight_quadric(q, magnitude(n_face), options.weighting);

        quadrics[face[(k + 1) % 3]] += q;
        quadrics[face[(k + 2) % 3]] += q;
      }
    }
  }
}

bool edge_topo_correctness(const F &faces, const E &edges) {
  for (const auto &face : faces) {
    for (order ord = 0; ord < 3; ++ord) {
      auto &vv = face.edge(ord)->vertices();
      assert(vv[0] < vv[1]);
      idx v_smaller = face[(ord + 1) % 3];
      idx v_larger = face[(ord + 2) % 3];
      if (v_smaller > v_larger) std::swap(v_smaller, v_larger);

      // edge vertices should match face corners
      assert(vv[0] == v_smaller && vv[1] == v_larger);
    }
  }

  return true;
}

void construct_edges(const V &vertices, F &faces, E &edges) {
  std::map<std::pair<idx, idx>, Edge> edge_set;

  // insert all edges into edge_set and find out if it is on boundary
  for (idx f = 0; f < faces.size(); ++f) {
    const auto &face = faces[f];
    for (order k = 0; k < 3; ++k) {
      // construct edge (v[i], v[j]);
      // edge local index will be k (= that of the 3rd vertex)
      const order i = (k + 1) % 3;
      const order j = (k + 2) % 3;
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
    faces[edge.face(0)].attach_edge(&edge, edge.ord_in_face(0));
    if (!edge.both_v_on_border())
      faces[edge.face(1)].attach_edge(&edge, edge.ord_in_face(1));
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
