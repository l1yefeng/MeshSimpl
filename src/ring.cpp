//
// Created by nickl on 5/13/19.
//

#include "ring.hpp"
#include <algorithm>     // for sort
#include <cassert>       // for assert
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
