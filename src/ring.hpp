//
// Created by nickl on 5/13/19.
//

#ifndef MESH_SIMPL_RING_HPP
#define MESH_SIMPL_RING_HPP

#include "neighbor.hpp"
#include "qem_heap.hpp"
#include "types.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

// 1-ring neighborhood of edge of interest, i.e., the union of 1-ring
// neighborhood of the endpoints of this edge.
class Ring {
 private:
  bool ccw_when_collect() const {
    return edge.ord_in_face(0) != (faces[edge.face(0)].v_order(v_del) + 1) % 3;
  }

 protected:
  static const size_t ESTIMATE_VALENCE = 8;

  // keep some references internally
  V &vertices;
  F &faces;
  std::vector<idx> v_del_twins, v_kept_twins;

  // in the center of the ring
  const Edge &edge;
  idx v_kept, v_del;

  // collections
  bool ccw;
  std::vector<Neighbor> v_del_neighbors, v_kept_neighbors;

  Ring(V &vertices, F &faces, const Edge &edge)
      : vertices(vertices),
        faces(faces),
        edge(edge),
        v_kept(edge[1 - edge.v_del_order()]),
        v_del(edge[edge.v_del_order()]),
        ccw(ccw_when_collect()) {}

  void reserve(size_t sz) {
    v_del_neighbors.reserve(sz);
    v_kept_neighbors.reserve(sz);
    v_del_twins.reserve(sz);
    v_kept_twins.reserve(sz);
  }

 public:
  virtual ~Ring() = default;

  virtual void collect() = 0;

  virtual bool check_env() = 0;

  bool check_topo() {
    return !sort_and_find_intersection(v_del_twins, v_kept_twins);
  }

  bool check_geom(double foldover_angle) const;

  bool check_quality(double aspect_ratio) const;

  virtual void collapse(Q &quadrics, QEMHeap &heap, bool fix_boundary) = 0;
};

class InteriorRing : public Ring {
 public:
  InteriorRing(V &vertices, F &faces, const Edge &edge)
      : Ring(vertices, faces, edge) {
    reserve(ESTIMATE_VALENCE);
  }

  void collect() override;

  bool check_env() override;

  void collapse(Q &quadrics, QEMHeap &heap, bool fix_boundary) override;
};

class BoundaryRing : public Ring {
 public:
  BoundaryRing(V &vertices, F &faces, const Edge &edge)
      : Ring(vertices, faces, edge) {
    reserve(ESTIMATE_VALENCE / 2);
  }

  void collect() override;

  bool check_env() override;

  void collapse(Q &quadrics, QEMHeap &heap, bool fix_boundary) override;
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_RING_HPP
