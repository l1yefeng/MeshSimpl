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
  const idx f = edge.face(0);

  Neighbor nb(f, edge.ordInF(0), ccw);
  while (true) {
    const Edge *nextEdge = nb.secondEdge(faces);

    if (nextEdge->onBoundary()) break;

    nb.rotate(faces);
    vDelNeighbors.push_back(nb);
  }

  nb = Neighbor(f, edge.ordInF(0), !ccw);
  while (true) {
    const Edge *nextEdge = nb.secondEdge(faces);

    if (nextEdge->onBoundary()) break;

    nb.rotate(faces);
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
    faces.setSide(it->f(), it->j(), edgeKept);
    dirtyEdge->replaceWing(f, it->f(), it->j());
  }

  faces.erase(f);

  // every face centered around deleted vertex
  for (const auto &nb : vDelNeighbors) {
    faces.setV(nb.f(), nb.center(), vKept);
    dirtyEdge = nb.secondEdge(faces);

    dirtyEdge->replaceEndpoint(vDel, vKept);

    assert(!dirtyEdge->neitherEndOnBoundary());
    updateEdge(dirtyEdge);
  }

  vertices.erase(vDel);

  updateEdge(dirtyEdge);
  for (auto nb : vKeptNeighbors) {
    dirtyEdge = nb.secondEdge(faces);
    updateEdge(dirtyEdge);
  }
}

}  // namespace Internal
}  // namespace MeshSimpl