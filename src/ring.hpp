//
// Created by nickl on 5/13/19.
//

#ifndef MESH_SIMPL_RING_HPP
#define MESH_SIMPL_RING_HPP

#include <cstddef>
#include <vector>
#include "edge.hpp"
#include "faces.hpp"
#include "neighbor.hpp"
#include "types.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

class QEMHeap;
class Vertices;

// 1-ring neighborhood of edge of interest, i.e., the union of 1-ring
// neighborhood of the endpoints of this edge.
class Ring {
 protected:
  static const size_t ESTIMATE_VALENCE = 8;

  // keep some references internally
  Vertices &vertices;
  Faces &faces;
  QEMHeap &heap;
  const SimplifyOptions &options;

  // in the center of the ring
  const Edge &edge;
  idx vKept, vDel;

  // collections
  std::vector<Neighbor> vDelNeighbors, vKeptNeighbors;

  Ring(Vertices &vertices, Faces &faces, QEMHeap &heap,
       const SimplifyOptions &options, const Edge &edge)
      : vertices(vertices),
        faces(faces),
        heap(heap),
        options(options),
        edge(edge),
        vKept(edge[1 - edge.delEndpointOrder()]),
        vDel(edge[edge.delEndpointOrder()]) {}

  void reserve(size_t sz) {
    vDelNeighbors.reserve(sz);
    vKeptNeighbors.reserve(sz);
  }

  virtual bool checkEnv() = 0;

  bool checkTopo();

  bool checkGeom(double foldOverAngle) const;

  bool checkQuality(double aspectRatio) const;

  void updateEdge(Edge *outdated);

 public:
  virtual ~Ring() = default;

  virtual void collect() = 0;

  bool check() {
    return checkEnv() && checkTopo() &&
           checkGeom(options.foldOverAngleThreshold) &&
           checkQuality(options.aspectRatioAtLeast);
  }

  virtual void collapse() = 0;
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_RING_HPP
