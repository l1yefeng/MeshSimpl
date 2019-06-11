//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_COLLAPSER_HPP
#define MESH_SIMPL_COLLAPSER_HPP

#include <array>
#include <cassert>
#include <initializer_list>
#include <vector>

#include "edge.hpp"
#include "faces.hpp"
#include "neighbor.hpp"
#include "proc.hpp"
#include "qemheap.hpp"
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Vertices;

class Collapser {
 private:
  Vertices& vertices;
  Faces& faces;
  QEMHeap& heap;
  Edge* target;
  const SimplifyOptions& options;

  std::array<std::vector<Neighbor>, 2> neighbors;
  int fRemoved;
  std::vector<Edge*> dirtyEdges;

  // Represent a pair of coincided edges. Although there is never a non-manifold
  // edge created during the whole process, the coincided edges will become
  // non-manifold in output if not handled beforehand thus the name.
  // NOTE: they are two different edges, with identical endpoints, certainly
  // different wings.
  struct NonManiInfo {
    idx vKept, vOther;
    std::array<Edge*, 2> edges;
    int status;
    NonManiInfo(idx vKept, idx vOther, Edge* e0, Edge* e1)
        : vKept(vKept), vOther(vOther), edges({e0, e1}), status(0) {}
  };
  std::vector<NonManiInfo> nonMani;

  void visitNonMani(idx vKept, idx vOther);

  void reset() {
    fRemoved = 0;
    target = nullptr;
    for (auto& n : neighbors) n.clear();
    dirtyEdges.clear();
    nonMani.clear();
  }

  int accept() {
    int temp = fRemoved;
    reset();
    return temp;
  }

  int reject() {
    heap.penalize(target);
    assert(fRemoved == 0);
    reset();
    return fRemoved;
  }

  void eraseF(idx f) {
    faces.erase(f);
    ++fRemoved;
  }

  // Store neighbors around endpoint(i) into neighbors[i], where i in {0, 1}
  void collect();

  // Find potential coincided edges (if collapse) and push to nonMani
  void findCoincideEdges(idx vKept);

  // Return true if no face will be flipped
  bool checkGeom() const {
    for (int i : {0, 1}) {
      for (const auto& nb : neighbors[i]) {
        if (isFaceFlipped(vertices, faces, nb.f(), nb.center(),
                          target->center(), options.foldOverAngleThreshold))
          return false;
      }
    }

    return true;
  }

  // Eliminate (not completely) coincided edges. This method will create a fork
  // for each endpoint of a pair of coincided edges and split the egdes by
  // assign the forked endpoints to one of them. As a result, one or more pair
  // of coincided edges (with common endpoint vKept) will be separated
  bool cleanup();

  void updateNonManiGroup(idx vKept, idx vFork);

 public:
  Collapser(Vertices& vertices, Faces& faces, QEMHeap& heap,
            const SimplifyOptions& options)
      : vertices(vertices),
        faces(faces),
        heap(heap),
        target(nullptr),
        options(options),
        neighbors(),
        fRemoved(0),
        nonMani() {}

  int collapse(Edge* edge);
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_COLLAPSER_HPP
