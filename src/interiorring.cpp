//
// Created by nickl on 5/31/19.
//

#include "interiorring.hpp"

namespace MeshSimpl {
namespace Internal {


void InteriorRing::collect() {
  const idx f0 = edge.face(0);
  const idx f1 = edge.face(1);

  // from f0 to f1 around vDel
  Neighbor nb(f0, edge.ordInF(0), ccw);
  while (true) {
    nb.rotate(faces);
    if (nb.f() == f1) break;
    vDelNeighbors.push_back(nb);
  }

  if (edge.neitherEndOnBoundary()) {
    // from f1 to f0 around vKept
    nb = Neighbor(f1, edge.ordInF(1), ccw);
    while (true) {
      nb.rotate(faces);
      if (nb.f() == f0) break;
      vKeptNeighbors.push_back(nb);
    }
  } else {
    // from f1 to boundary around vKept
    nb = Neighbor(f1, edge.ordInF(1), ccw);
    while (!nb.secondEdge(faces)->onBoundary()) {
      nb.rotate(faces);
      vKeptNeighbors.push_back(nb);
    }

    // from f0 to boundary around vKept (switch direction)
    nb = Neighbor(f0, edge.ordInF(0), !ccw);
    while (!nb.secondEdge(faces)->onBoundary()) {
      nb.rotate(faces);
      vKeptNeighbors.push_back(nb);
    }
  }
}
bool InteriorRing::checkEnv() {
  // special case: there are 2 faces, 3 vertices in current component
  // this happens if input contains such component because
  // this edge collapse implementation avoid generating them
  if (vDelNeighbors.empty()) {
    assert(faces[edge.face(0)][edge.ordInF(0)] ==
           faces[edge.face(1)][edge.ordInF(1)]);
    return false;
  }

  // special case: avoid collapsing a tetrahedron into folded faces
  if (vDelNeighbors.size() == 1 && vKeptNeighbors.size() == 1) return false;

  // fine target to collapse
  return true;
}
void InteriorRing::collapse() {
  const idx f0 = edge.face(0);
  const idx f1 = edge.face(1);

  // update vertex position and quadric
  vertices.setPosition(vKept, edge.center());
  vertices.setQ(vKept, edge.q());

  // two kept edges in two deleted faces
  Edge *edgeKept0 = faces.edgeAcrossFrom(f0, vDel);
  Edge *edgeKept1 = faces.edgeAcrossFrom(f1, vDel);

  auto it = vDelNeighbors.begin();
  Edge *dirtyEdge = edgeKept0;

  // first face to process: deleted face 0
  faces.setV(it->f(), it->center(), vKept);
  faces.side(it->f(), it->j())->erase();

  faces.setSide(it->f(), it->j(), edgeKept0);
  dirtyEdge->replaceWing(f0, it->f(), it->j());
  faces.erase(f0);

  // every face centered around deleted vertex
  for (++it; it != vDelNeighbors.end(); ++it) {
    faces.setV(it->f(), it->center(), vKept);
    dirtyEdge = faces.side(it->f(), it->j());

    dirtyEdge->replaceEndpoint(vDel, vKept);

    updateEdge(dirtyEdge);
  }

  // the deleted face 1
  --it;
  faces.side(it->f(), it->i())->erase();
  faces.setSide(it->f(), it->i(), edgeKept1);
  dirtyEdge = edgeKept1;
  dirtyEdge->replaceWing(f1, it->f(), it->i());
  faces.erase(f1);

  vertices.erase(vDel);

  // every edge centered around the kept vertex
  for (auto nb : vKeptNeighbors) {
    dirtyEdge = nb.firstEdge(faces);
    updateEdge(dirtyEdge);
  }

  updateEdge(dirtyEdge);
}

}
}  // namespace MeshSimpl