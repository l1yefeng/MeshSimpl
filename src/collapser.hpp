//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_COLLAPSER_HPP
#define MESH_SIMPL_COLLAPSER_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "edge.hpp"
#include "faces.hpp"
#include "neighbor.hpp"
#include "proc.hpp"
#include "qemheap.hpp"
#include "types.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

class Collapser {
 private:
  Vertices& vertices;
  Faces& faces;
  QEMHeap& heap;
  Edge* target;
  const SimplifyOptions& options;

  std::array<std::vector<Neighbor>, 2> neighbors;
  int vRemoved;
  int fRemoved;

  typedef std::map<idx, std::array<Edge*, 2>> NonManiInfo;
  typedef std::pair<int, int> ReturnType;

  ReturnType accept() { return {vRemoved, fRemoved}; }

  ReturnType reject() {
    heap.penalize(target);
    assert(vRemoved == 0);
    assert(fRemoved == 0);
    return {vRemoved, fRemoved};
  }

  void eraseV(idx v) {
    vertices.erase(v);
    ++vRemoved;
  }

  idx duplicateV(idx src) {
    --vRemoved;
    return vertices.duplicate(src);
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

  void cleanup(std::map<idx, NonManiInfo>& nonManiGroup);

 public:
  Collapser(Vertices& vertices, Faces& faces, QEMHeap& heap, Edge* target,
            const SimplifyOptions& options)
      : vertices(vertices),
        faces(faces),
        heap(heap),
        target(target),
        options(options),
        neighbors(),
        vRemoved(0),
        fRemoved(0) {}

  ReturnType collapse();
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_COLLAPSER_HPP
