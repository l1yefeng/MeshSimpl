//
// Created by nickl on 5/31/19.
//

#include <cassert>
#include <vector>

#include "boundaryring.hpp"
#include "edge.hpp"
#include "faces.hpp"
#include "neighbor.hpp"
#include "types.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

void BoundaryRing::collect() {
  Neighbor nb(&edge, 0, vDel, faces);
  while (true) {
    const Edge *nextEdge = nb.secondEdge();

    if (nextEdge->onBoundary()) break;

    nb.rotate();
    vDelNeighbors.push_back(nb);
  }

  nb.replace(&edge, 0, vKept);
  while (true) {
    const Edge *nextEdge = nb.secondEdge();

    if (nextEdge->onBoundary()) break;

    nb.rotate();
    vKeptNeighbors.push_back(nb);
  }
}

bool BoundaryRing::checkEnv() {
  // special case: there is ONE triangle in current component
  if (vKeptNeighbors.empty()) return false;

  return true;
}

void BoundaryRing::collapse() {
  const idx f = edge.face(0);

  // update vertex position and quadric
  vertices.setPosition(vKept, edge.center());
  vertices.setQ(vKept, edge.q());

  Edge *edgeKept = faces.edgeAcrossFrom(f, vDel);
  Edge *dirtyEdge = edgeKept;

  // it is possible that vDelNeighbors is empty
  faces.edgeAcrossFrom(f, vKept)->erase();
  if (vDelNeighbors.empty()) {
    dirtyEdge->dropWing(f);
  } else {
    auto it = vDelNeighbors.begin();
    faces.setSide(it->f(), it->secondVOrd(), edgeKept);
    dirtyEdge->replaceWing(f, it->f(), it->secondVOrd());
  }

  faces.erase(f);

  // every face centered around deleted vertex
  for (const auto &nb : vDelNeighbors) {
    faces.setV(nb.f(), nb.center(), vKept);
    dirtyEdge = nb.secondEdge();

    dirtyEdge->replaceEndpoint(vDel, vKept);

    assert(!dirtyEdge->neitherEndOnBoundary());
    updateEdge(dirtyEdge);
  }

  vertices.erase(vDel);

  updateEdge(dirtyEdge);
  for (auto nb : vKeptNeighbors) {
    dirtyEdge = nb.secondEdge();
    updateEdge(dirtyEdge);
  }
}

}  // namespace Internal
}  // namespace MeshSimpl