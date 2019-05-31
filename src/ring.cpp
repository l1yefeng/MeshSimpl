//
// Created by nickl on 5/13/19.
//

#include "ring.hpp"
#include <algorithm>  // for sort
#include <cassert>    // for assert
#include "boundaryring.hpp"
#include "interiorring.hpp"
#include "proc.hpp"      // for updateEdge, isFaceElongated
#include "qemheap.hpp"   // for QEMHeap
#include "vertices.hpp"  // for Vertices

namespace MeshSimpl {
namespace Internal {

bool Ring::checkTopo() {
  std::vector<const Neighbor *> starDel, starKept;
  starDel.reserve(vDelNeighbors.size());
  starKept.reserve(vKeptNeighbors.size());

  for (const auto &nb : vDelNeighbors) starDel.push_back(&nb);
  for (const auto &nb : vKeptNeighbors) starKept.push_back(&nb);

  const auto cmp = [&](const Neighbor *nb0, const Neighbor *nb1) -> bool {
    return nb0->secondV(faces) < nb1->secondV(faces);
  };
  std::sort(starDel.begin(), starDel.end(), cmp);
  std::sort(starKept.begin(), starKept.end(), cmp);
  for (auto itd = starDel.begin(), itk = starKept.begin();
       itd != starDel.end() && itk != starKept.end();) {
    if (cmp(*itd, *itk)) {
      ++itd;
    } else if (cmp(*itk, *itd)) {
      ++itk;
    } else {
      return false;
    }
  }

  return true;
}

bool Ring::checkGeom(double foldOverAngle) const {
  for (const auto &nb : vDelNeighbors) {
    if (isFaceFolded(vertices, faces, nb.f(), nb.center(), edge.center(),
                     foldOverAngle))
      return false;
  }
  for (const auto &nb : vKeptNeighbors) {
    if (isFaceFolded(vertices, faces, nb.f(), nb.center(), edge.center(),
                     foldOverAngle))
      return false;
  }
  return true;
}

bool Ring::checkQuality(double aspectRatio) const {
  for (const auto &nb : vDelNeighbors) {
    if (isFaceElongated(vertices[vDel], vertices[nb.firstV(faces)],
                        vertices[nb.secondV(faces)], aspectRatio))
      return false;
  }
  for (const auto &nb : vDelNeighbors) {
    if (isFaceElongated(vertices[vKept], vertices[nb.firstV(faces)],
                        vertices[nb.secondV(faces)], aspectRatio))
      return false;
  }
  return true;
}

void Ring::updateEdge(Edge *outdated) {
  const double errorPrev = outdated->error();
  outdated->planCollapse(options.fixBoundary);
  heap.fix(outdated, errorPrev);
}

}  // namespace Internal
}  // namespace MeshSimpl
