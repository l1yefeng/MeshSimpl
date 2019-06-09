//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_COLLAPSER_HPP
#define MESH_SIMPL_COLLAPSER_HPP

#include <array>
#include <cassert>
#include <initializer_list>
#include <map>
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

  typedef std::map<idx, std::array<Edge*, 2>> NonManiInfo;
  typedef std::map<idx, NonManiInfo> NonManiGroup;

  NonManiGroup nonManiGroup;

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

  void findCoincideEdges(NonManiInfo& nonMani);

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

  void cleanup();

  void updateNonManiGroup(const std::map<idx, int>& visited,
                          NonManiGroup::iterator current, idx vFork);

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
        nonManiGroup() {}

  int collapse();
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_COLLAPSER_HPP
