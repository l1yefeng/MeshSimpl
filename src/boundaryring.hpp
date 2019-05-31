//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_BOUNDARYRING_HPP
#define MESH_SIMPL_BOUNDARYRING_HPP

#include "ring.hpp"
namespace MeshSimpl {
namespace Internal {

class BoundaryRing : public Ring {
 public:
  BoundaryRing(Vertices &vertices, Faces &faces, QEMHeap &heap,
               const SimplifyOptions &options, const Edge &edge)
      : Ring(vertices, faces, heap, options, edge) {
    reserve(ESTIMATE_VALENCE / 2);
  }

  void collect() override;

  bool checkEnv() override;

  void collapse() override;
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_BOUNDARYRING_HPP
