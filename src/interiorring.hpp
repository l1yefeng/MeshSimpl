//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_INTERIORRING_HPP
#define MESH_SIMPL_INTERIORRING_HPP

#include "ring.hpp"

namespace MeshSimpl {

struct SimplifyOptions;

namespace Internal {

class Edge;
class Faces;
class Vertices;
class QEMHeap;

class InteriorRing : public Ring {
 public:
  InteriorRing(Vertices &vertices, Faces &faces, QEMHeap &heap,
               const SimplifyOptions &options, const Edge &edge)
      : Ring(vertices, faces, heap, options, edge) {
    reserve(ESTIMATE_VALENCE);
  }

  void collect() override;

  bool checkEnv() override;

  void collapse() override;
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_INTERIORRING_HPP
