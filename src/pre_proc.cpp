//
// Created by nickl on 1/8/19.
//

#include "pre_proc.hpp"
#include <algorithm>
#include <limits>
#include "edge.hpp"
#include "marker.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

void weight_quadric(Quadric& quadric, double face_area, WEIGHTING strategy) {
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

void compute_quadrics(const V& vertices, Connectivity& conn, Q& quadrics,
                      const SimplifyOptions& options) {
  // quadrics are initialized with all zeros
  quadrics.resize(vertices.size());

  for (const auto& face : conn.indices) {
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

    for (const auto& v : face) quadrics[v] += q;
  }

  if (!options.fix_boundary) {
    for (idx f = 0; f < conn.indices.size(); ++f) {
      order k;
      for (k = 0; k < 3; ++k)
        if (conn.edge_of_face(f, k).on_boundary()) break;

      if (k == 3) continue;

      const auto& v = conn.indices[f];
      const std::array<vec3d, 3> e{vertices[v[2]] - vertices[v[1]],
                                   vertices[v[0]] - vertices[v[2]],
                                   vertices[v[1]] - vertices[v[0]]};
      const vec3d n_face = cross(e[0], e[1]);

      for (k = 0; k < 3; ++k) {
        if (!conn.edge_of_face(f, k).on_boundary()) continue;

        vec3d normal = cross(n_face, e[k]);
        const double normal_mag = magnitude(normal);
        if (normal_mag != 0)
          normal /= normal_mag;
        else
          continue;

        const double d = -dot(normal, vertices[v[(k + 1) % 3]]);
        Quadric q = make_quadric(normal, d);
        q *= options.border_constraint;
        weight_quadric(q, magnitude(n_face), options.weighting);

        quadrics[v[(k + 1) % 3]] += q;
        quadrics[v[(k + 2) % 3]] += q;
      }
    }
  }
}

bool edge_topo_correctness(const Connectivity& conn) {
  for (idx f = 0; f < conn.indices.size(); ++f) {
    const auto& f2e = conn.face2edge[f];
    for (order i = 0; i < 3; ++i) {
      auto& vv = conn.edges[f2e[i]].vertices();
      assert(vv[0] < vv[1]);
      idx v_smaller = conn.indices[f][(i + 1) % 3];
      idx v_larger = conn.indices[f][(i + 2) % 3];
      if (v_smaller > v_larger) std::swap(v_smaller, v_larger);

      // edge vertices should match face corners
      assert(vv[0] == v_smaller && vv[1] == v_larger);
    }
  }

  return true;
}

void construct_edges(const V& vertices, Internal::Connectivity& conn) {
  std::map<std::pair<idx, idx>, Edge> edge_set;

  // insert all edges into edge_set and find out if it is on boundary
  for (idx f = 0; f < conn.indices.size(); ++f) {
    const auto& face = conn.indices[f];
    for (order k = 0; k < 3; ++k) {
      // construct edge (v[i], v[j]);
      // edge local index will be k (= that of the 3rd vertex)
      const order i = (k + 1) % 3;
      const order j = (k + 2) % 3;
      idx v0 = face[i];
      idx v1 = face[j];
      if (v0 > v1) std::swap(v0, v1);
      Edge edge(v0, v1);
      edge.attach_1st_face(f, k);
      auto it_and_inserted = edge_set.emplace(std::make_pair(v0, v1), edge);
      auto it = it_and_inserted.first;

      if (!it_and_inserted.second) {
        if (!it->second.attach_2nd_face(f, k)) {
          throw std::invalid_argument(
              "ERROR::INPUT_MESH: detected non-manifold edge");
        }
      }
    }
  }

  // populate edges vector from map
  conn.edges.reserve(edge_set.size());
  for (const auto& elem : edge_set) conn.edges.emplace_back(elem.second);

  conn.face2edge.resize(conn.indices.size());
  std::vector<bool> vertex_on_boundary(vertices.size(), false);

  for (idx e = 0; e < conn.edges.size(); ++e) {
    const auto& edge = conn.edges[e];
    // identify boundary vertices
    if (edge.both_v_on_border()) {
      vertex_on_boundary[edge[0]] = true;
      vertex_on_boundary[edge[1]] = true;
    }

    // populate face2edge references
    conn.face2edge[edge.face(0)][edge.ord_in_face(0)] = e;
    if (!edge.both_v_on_border())
      conn.face2edge[edge.face(1)][edge.ord_in_face(1)] = e;
  }

  // non-boundary edges may have one vertex on boundary, find them in this loop
  for (auto& edge : conn.edges) {
    edge.set_v_on_border(vertex_on_boundary[edge[0]],
                         vertex_on_boundary[edge[1]]);
  }

  assert(edge_topo_correctness(conn));
}

}  // namespace Internal
}  // namespace MeshSimpl
