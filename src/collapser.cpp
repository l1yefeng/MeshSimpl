//
// Created by nickl on 5/31/19.
//

#include <algorithm>
#include <tuple>
#include <utility>

#include "collapser.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

void Collapser::collect() {
  for (order i : {0, 1}) {
    idx v = target->endpoint(i);
    if (!vertices.isBoundary(v)) {
      // traverse around like a fan
      Neighbor nb(target, i, v, faces);
      for (nb.rotate(); nb.f() != target->face(1 - i); nb.rotate()) {
        neighbors[i].push_back(nb);
      }
    } else {
      // traverse with two rows stop at boundary
      for (order column : {0, 1}) {
        Neighbor nb(target, column, v, faces);
        while (!nb.secondEdge()->onBoundary()) {
          nb.rotate();
          neighbors[i].push_back(nb);
        }

        // if true, target only has face[0]
        if (target->onBoundary()) break;
      }
    }
  }
}

void Collapser::findCoincideEdges(NonManiInfo& nonMani) {
  // comparision used in sorting and finding intersection
  const auto cmp = [](const Neighbor* nb0, const Neighbor* nb1) -> bool {
    return nb0->secondV() < nb1->secondV();
  };

  std::array<std::vector<const Neighbor*>, 2> pointers;
  for (int i : {0, 1}) {
    pointers[i].reserve(neighbors[i].size());
    for (const auto& nb : neighbors[i]) pointers[i].push_back(&nb);
    std::sort(pointers[i].begin(), pointers[i].end(), cmp);
  }

  // finding intersection
  for (auto it0 = pointers[0].begin(), it1 = pointers[1].begin();
       it0 != pointers[0].end() && it1 != pointers[1].end();) {
    if (cmp(*it0, *it1)) {
      ++it0;
    } else if (cmp(*it1, *it0)) {
      ++it1;
    } else {
      std::array<Edge*, 2> coincided = {(*it0)->secondEdge(),
                                        (*it1)->secondEdge()};
      nonMani.emplace((*it0)->secondV(), coincided);

      ++it0;
      ++it1;
    }
  }
}

void Collapser::cleanup(std::map<idx, NonManiInfo>& nonManiGroup) {
  assert(!nonManiGroup.empty());
  // the original vKept or a fork
  idx vKept = nonManiGroup.begin()->first;
  NonManiInfo& nonMani = nonManiGroup.begin()->second;
  assert(!nonMani.empty());

  // use the first non-manifold edge to separate this mess
  std::vector<std::tuple<Edge*, idx, idx>> edgesReplaceEnd;
  std::vector<std::tuple<idx, order, idx>> facesSetV;

  auto& nm = *nonMani.begin();
  // e0 will be attached to the e0->face(0) now and the face connects to it
  // on e1 e0 will keep current endpoints e1 will be attached to the
  // e0->face(1) now and the face connects to it on e1 e1 will have
  // endpoints vKept replaced with vDel, vOther with vOtherFork
  Edge* e0 = nm.second[0];
  Edge* e1 = nm.second[1];
  assert(e0->endpoints() == e1->endpoints());

  idx vKeptFork = duplicateV(vKept);
  idx vOther = nm.first;
  idx vOtherFork = duplicateV(vOther);
  vertices.setBoundary(vKeptFork, false);

  idx fExch0;

  // used to check whether another non-manifold edge will be handled
  // while we fork this non-manifold edge (visited once meaning yeah)
  std::map<idx, int> visited;
  for (auto& _nm : nonMani) {
    visited.emplace(_nm.first, 0);
  }

  // select a direction (1, but 0 will work as well), circle around
  // until met the coincided edge and every visited face will be separated
  // to the forked edge to turn the non-manifold into 2-manifold
  order traverseOrd = 0;
  Neighbor nb(e0, traverseOrd, vKept, faces);
  assert(faces.v(nb.f(), nb.center()) == vKept);
  while (true) {
    edgesReplaceEnd.emplace_back(nb.secondEdge(), vKept, vKeptFork);
    facesSetV.emplace_back(nb.f(), nb.center(), vKeptFork);

    auto itVisited = visited.find(nb.secondV());
    if (itVisited != visited.end()) {
      itVisited->second += 1;
    }

    if (nb.secondEdge() == e1) {
      fExch0 = nb.f();
      break;
    }
    if (nb.secondEdge()->onBoundary()) {
      // switch direction in order to separate edge in one pass
      edgesReplaceEnd.clear();
      facesSetV.clear();
      for (auto& vis : visited) vis.second = 0;

      traverseOrd = 1 - traverseOrd;
      nb.replace(e0, traverseOrd, vKept);
      continue;
    }

    nb.rotate();
    assert(nb.secondEdge() != e0);
  }

  nb.replace(e0, traverseOrd, vOther);
  assert(e0->ordInF(traverseOrd) != INVALID);
  assert(faces.v(nb.f(), nb.center()) == vOther);
  bool hitBoundary = false;
  while (true) {
    edgesReplaceEnd.emplace_back(nb.secondEdge(), vOther, vOtherFork);
    facesSetV.emplace_back(nb.f(), nb.center(), vOtherFork);

    if (hitBoundary) assert(nb.secondEdge() != e1);

    if (nb.secondEdge() == e1) {
      assert(fExch0 == nb.f());
      assert(!hitBoundary);
      break;
    }
    if (nb.secondEdge()->onBoundary()) {
      if (!hitBoundary) {
        nb.replace(e1, e1->wingOrder(fExch0), vOther);
        hitBoundary = true;
        continue;
      } else {
        edgesReplaceEnd.emplace_back(e1, vOther, vOtherFork);
        break;
      }
    }

    nb.rotate();
    assert(nb.secondEdge() != e0);
  }
  vertices.setBoundary(vOtherFork, hitBoundary);

  // exchange wing between e0 and e1. when completed, e0 will have
  // two wings remain on surface while e1 will have two wings
  // detached (endpoints are forked)
  order fExch1Ord = 1 - e1->wingOrder(fExch0);
  faces.setSide(e0->face(traverseOrd), e0->ordInF(traverseOrd), e1);
  if (!e1->onBoundary()) {
    faces.setSide(e1->face(fExch1Ord), e1->ordInF(fExch1Ord), e0);
    idx fExch1 = e1->face(fExch1Ord);
    order e1OrdInFExch = e1->ordInF(fExch1Ord);
    e1->setWing(fExch1Ord, e0->face(traverseOrd), e0->ordInF(traverseOrd));
    e0->setWing(traverseOrd, fExch1, e1OrdInFExch);
  } else {
    assert(fExch1Ord == 1);
    e1->setWing(fExch1Ord, e0->face(traverseOrd), e0->ordInF(traverseOrd));
    bool edgeValid = e0->dropWing(e0->face(traverseOrd));
    if (!edgeValid) {
      e0->erase();
      eraseV(e0->endpoint(0));
      eraseV(e0->endpoint(1));
    }
  }

  for (auto& ere : edgesReplaceEnd)
    std::get<0>(ere)->replaceEndpoint(std::get<1>(ere), std::get<2>(ere));
  for (auto& fsv : facesSetV)
    faces.setV(std::get<0>(fsv), std::get<1>(fsv), std::get<2>(fsv));

  for (auto it = nonMani.begin(), last = nonMani.end(); it != last;) {
    switch (visited[it->first]) {
      case 1:
        it = nonMani.erase(it);
        break;
      case 0:
        ++it;
        break;
      case 2: {
        auto insertRes = nonManiGroup.emplace(vKeptFork, NonManiInfo{});
        insertRes.first->second.emplace(it->first, it->second);
        it = nonMani.erase(it);
        break;
      }
      default:
        assert(false);
    }
  }

  if (nonMani.empty()) {
    nonManiGroup.erase(vKept);
  }
}

Collapser::ReturnType Collapser::collapse() {
  collect();

  order delOrd = vertices.isBoundary(target->endpoint(0)) ? 1 : 0;

  bool neck = target->bothEndsOnBoundary() && !target->onBoundary();

  // check cause of topo change
  idx vDel = target->endpoint(delOrd);
  idx vKept = target->endpoint(1 - delOrd);
  NonManiInfo nonMani;
  findCoincideEdges(nonMani);

  // reject now, before any modification that changes topology is applied
  if (!options.topologyModifiable) {
    // collapse will create non-manifold edges
    if (!nonMani.empty()) {
      return reject();
    }

    // reject if it will erase/separate the component
    if (neck) {
      std::array<bool, 2> single{false, false};
      for (order i : {0, 1}) {
        idx f = target->face(i);
        order ord = target->ordInF(i);
        single[i] = faces.side(f, next(ord))->onBoundary() &&
                    faces.side(f, prev(ord))->onBoundary();
      }

      if (single[0] == single[1]) {
        return reject();
      }
    }

    // reject if what remains is a tetrahedron
    if (target->bothEndsOnBoundary() && neighbors[0].size() == 1 &&
        neighbors[1].size() == 1) {
      return reject();
    }

    // reject if what remains are two faces folded back-to-back
    if (!vertices.isBoundary(vDel) && neighbors[delOrd].empty()) {
      return reject();
    }
  }

  // reject if topology preserves but some face will be flipped
  if (nonMani.empty()) {
    if (!checkGeom()) {
      return reject();
    }
  }

  // special case: two faces folded (#f=2, #v=3)
  if (!vertices.isBoundary(vDel) && neighbors[delOrd].empty()) {
    idx f0 = target->face(0);
    idx f1 = target->face(1);
    for (order ord : {0, 1, 2}) {
      Edge* edge = faces.side(f0, ord);
      assert(edge->face(0) + edge->face(1) == f0 + f1);

      edge->erase();
      eraseV(faces.v(f0, ord));
    }
    eraseF(f0);
    eraseF(f1);
    return accept();
  }

  // maps from center vertex to a group of coincide edges
  std::map<idx, NonManiInfo> nonManiGroup{};

  // there is topo change or not, collapse the target now. cleanup afterwords
  // update vertex data
  vertices.setPosition(vKept, target->center());
  vertices.setQ(vKept, target->q());

  // replace face corner
  for (auto& nb : neighbors[delOrd]) {
    faces.setV(nb.f(), nb.center(), vKept);
  }

  std::array<std::vector<Edge*>, 2> dirtyEdges;
  // collect edges who need update around endpoint 0 and 1
  for (int i : {0, 1}) {
    if (!vertices.isBoundary(target->endpoint(i))) {
      auto it = neighbors[i].begin();
      dirtyEdges[i].push_back(it->firstEdge());
      for (; it != neighbors[i].end(); ++it) {
        dirtyEdges[i].push_back(it->secondEdge());
      }
    } else {
      for (int column : {0, 1}) {
        dirtyEdges[i].push_back(faces.edgeAcrossFrom(target->face(column),
                                                     target->endpoint(1 - i)));
        if (target->onBoundary()) break;
      }
      for (auto& nb : neighbors[i]) {
        dirtyEdges[i].push_back(nb.secondEdge());
      }
    }
  }

  // replace edge endpoint
  for (auto& edge : dirtyEdges[delOrd]) {
    edge->replaceEndpoint(vDel, vKept);
  }

  // update error of edges
  dirtyEdges[0].insert(dirtyEdges[0].end(), dirtyEdges[1].begin(),
                       dirtyEdges[1].end());
  for (auto& dirty : dirtyEdges[0]) {
    if (options.fixBoundary && dirty->bothEndsOnBoundary()) {
      heap.remove(dirty);
      continue;
    }
    double errorPrev = dirty->error();
    dirty->planCollapse(options.fixBoundary);
    heap.fix(dirty, errorPrev);
  }

  // take away face 0 and 1
  std::array<bool, 2> edgeValid = {true, true};
  std::array<Edge*, 2> edgeKept = {nullptr, nullptr};
  for (int i : {0, 1}) {
    idx fDel = target->face(i);
    edgeKept[i] = faces.edgeAcrossFrom(fDel, vDel);
    Edge* edgeDel = faces.edgeAcrossFrom(fDel, vKept);

    if (!edgeDel->onBoundary()) {
      order wingOrd = 1 - edgeDel->wingOrder(fDel);
      idx fDirty = edgeDel->face(wingOrd);
      idx edgeDelOrdInFDirty = edgeDel->ordInF(wingOrd);
      edgeKept[i]->setWing(edgeKept[i]->wingOrder(fDel), fDirty,
                           edgeDelOrdInFDirty);
      faces.setSide(fDirty, edgeDelOrdInFDirty, edgeKept[i]);
    } else {
      edgeValid[i] = edgeKept[i]->dropWing(fDel);
      if (!edgeValid[i]) {
        eraseV(faces.v(fDel, target->ordInF(i)));
        edgeKept[i]->erase();
      }
    }
    edgeDel->erase();

    eraseF(fDel);

    if (target->onBoundary()) break;
  }

  target->erase();
  eraseV(vDel);

  // special case: one face or two faces with one common edge
  if (!edgeValid[0] && (!edgeKept[1] || !edgeValid[1])) {
    eraseV(vKept);
    assert(vRemoved == (edgeKept[1] ? 4 : 3));
    return accept();
  }

  // special case: component is separated
  // FIXME: ERROR when nonMani has one of the coincided edges updated
  if (neck && edgeValid[0] && edgeValid[1]) {
    Edge* seed = edgeKept[1];
    idx vFork = duplicateV(vKept);
    assert(nonMani.find(seed->endpoint(0) + seed->endpoint(1) - vFork) ==
           nonMani.end());
    std::vector<Neighbor> dirtyNeighbors;
    for (int column : {0, 1}) {
      Neighbor nb(seed, column, vKept, faces);
      while (true) {
        dirtyNeighbors.push_back(nb);
        NonManiInfo::iterator it;
        if ((it = nonMani.find(nb.secondV())) != nonMani.end()) {
          auto res = nonManiGroup.emplace(vFork, NonManiInfo{});
          res.first->second.emplace(it->first, it->second);
          nonMani.erase(it);
        }
        if (nb.secondEdge()->onBoundary()) break;
        nb.rotate();
      }
      if (seed->onBoundary()) break;
    }

    seed->replaceEndpoint(vKept, vFork);
    for (auto& nb : dirtyNeighbors) {
      faces.setV(nb.f(), nb.center(), vFork);
      nb.secondEdge()->replaceEndpoint(vKept, vFork);
    }

    assert(vRemoved == 0);
  }

  if (!nonMani.empty()) {
    nonManiGroup.emplace(vKept, nonMani);
  }

  while (!nonManiGroup.empty()) {
    cleanup(nonManiGroup);
  }

  return accept();
}

}  // namespace Internal
}