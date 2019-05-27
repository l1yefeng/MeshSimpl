//
// Created by nickl on 5/13/19.
//

#ifndef MESH_SIMPL_RING_HPP
#define MESH_SIMPL_RING_HPP

#include <cstddef>       // for size_t
#include <vector>        // for vector
#include "edge.hpp"      // for Edge
#include "faces.hpp"     // for Faces
#include "neighbor.hpp"  // for Neighbor
#include "types.hpp"     // for SimplifyOptions, idx
#include "util.hpp"      // for next

namespace MeshSimpl {
namespace Internal {

class QEMHeap;
class Vertices;

// 1-ring neighborhood of edge of interest, i.e., the union of 1-ring
// neighborhood of the endpoints of this edge.
class Ring {
 private:
  bool ccwWhenCollect() const {
    return edge.ordInF(0) != next(faces.orderOf(edge.face(0), vDel));
  }

 protected:
  static const size_t ESTIMATE_VALENCE = 8;

  // keep some references internally
  Vertices &vertices;
  Faces &faces;

  // in the center of the ring
  const Edge &edge;
  idx vKept, vDel;

  // collections
  bool ccw;
  std::vector<Neighbor> vDelNeighbors, vKeptNeighbors;

  Ring(Vertices &vertices, Faces &faces, const Edge &edge)
      : vertices(vertices),
        faces(faces),
        edge(edge),
        vKept(edge[1 - edge.delEndpointOrder()]),
        vDel(edge[edge.delEndpointOrder()]),
        ccw(ccwWhenCollect()) {}

  void reserve(size_t sz) {
    vDelNeighbors.reserve(sz);
    vKeptNeighbors.reserve(sz);
  }

  virtual bool checkEnv() = 0;

  bool checkTopo();

  bool checkGeom(double foldOverAngle) const;

  bool checkQuality(double aspectRatio) const;

 public:
  virtual ~Ring() = default;

  virtual void collect() = 0;

  bool check(const SimplifyOptions &options) {
    return checkEnv() && checkTopo() &&
           checkGeom(options.foldOverAngleThreshold) &&
           checkQuality(options.aspectRatioAtLeast);
  }

  virtual void collapse(QEMHeap &heap, bool fixBoundary) = 0;
};

class InteriorRing : public Ring {
 public:
  InteriorRing(Vertices &vertices, Faces &faces, const Edge &edge)
      : Ring(vertices, faces, edge) {
    reserve(ESTIMATE_VALENCE);
  }

  void collect() override;

  bool checkEnv() override;

  void collapse(QEMHeap &heap, bool fixBoundary) override;
};

class BoundaryRing : public Ring {
 public:
  BoundaryRing(Vertices &vertices, Faces &faces, const Edge &edge)
      : Ring(vertices, faces, edge) {
    reserve(ESTIMATE_VALENCE / 2);
  }

  void collect() override;

  bool checkEnv() override;

  void collapse(QEMHeap &heap, bool fixBoundary) override;
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_RING_HPP
