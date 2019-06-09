//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_COLLAPSER_HPP
#define MESH_SIMPL_COLLAPSER_HPP

#include <algorithm>
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

  struct NonManiInfo {
    idx vKept, vOther;
    std::array<Edge*, 2> edges;
    int status;
    NonManiInfo(idx vKept, idx vOther, Edge* e0, Edge* e1)
        : vKept(vKept), vOther(vOther), edges({e0, e1}), status(0) {}
  };
  std::vector<NonManiInfo> nonMani;

  void visitNonMani(idx vKept, idx vOther) {
    for (auto& nm : nonMani) {
      if (nm.status >= 0 && nm.vKept == vKept && nm.vOther == vOther) {
        ++nm.status;
        assert(nm.status == 1 || nm.status == 2);
        return;
      }
    }
  }

  std::vector<NonManiInfo>::iterator frontNM() {
    return std::find_if(nonMani.begin(), nonMani.end(),
                        [](const NonManiInfo& nm) { return nm.status >= 0; });
  }

  int accept() { return fRemoved; }

  int reject() {
    heap.penalize(target);
    assert(fRemoved == 0);
    return fRemoved;
  }

  void eraseF(idx f) {
    faces.erase(f);
    ++fRemoved;
  }

  // Store neighbors around endpoint(i) into neighbors[i], where i in {0, 1}
  void collect();

  void findCoincideEdges(idx vKept);

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

  bool cleanup();

  void updateNonManiGroup(idx vKept, idx vFork);

 public:
  Collapser(Vertices& vertices, Faces& faces, QEMHeap& heap, Edge* target,
            const SimplifyOptions& options)
      : vertices(vertices),
        faces(faces),
        heap(heap),
        target(target),
        options(options),
        neighbors(),
        fRemoved(0),
        nonMani() {}

  int collapse();
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_COLLAPSER_HPP
