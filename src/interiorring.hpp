//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_INTERIORRING_HPP
#define MESH_SIMPL_INTERIORRING_HPP

#include "ring.hpp"
namespace MeshSimpl {
namespace Internal {

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



}
}  // namespace MeshSimpl

#endif //MESH_SIMPL_INTERIORRING_HPP
