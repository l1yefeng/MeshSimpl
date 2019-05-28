//
// Created by nickl on 1/8/19.
//

#include "qemheap.hpp"
#include <cassert>   // for assert
#include <cmath>     // for isnan
#include <memory>    // for allocator_traits<>::value_type
#include <utility>   // for swap
#include "edge.hpp"  // for Edge

namespace MeshSimpl {
namespace Internal {

QEMHeap::QEMHeap(Edges &edges )
    : keys(edges.size() + 1), edges(edges), handles(edges.size()), n(0) {
}

void QEMHeap::pop() {
  exchange(1, n--);
  sink(1);
  keys.resize(n + 1);
}

void QEMHeap::fix(const Edge *ptr, double errorPrev) {
  auto e = static_cast<idx>(ptr - edges.data());
  size_t k = handles[e];
  if (ptr->error() > errorPrev)
    sink(k);
  else
    swim(k);
}

void QEMHeap::penalize(Edge *edge) {
  edge->setErrorInfty();
  sink(handles[edge - edges.data()]);
}

void QEMHeap::erase(const Edge *ptr) {
  idx e = ptr - edges.data();

  const double errorPrev = edges[e].error();
  size_t k = handles[e];

  if (k < n + 1) {
    exchange(k, n--);
    if (edges[keys[k]].error() > errorPrev)
      sink(k);
    else
      swim(k);
    keys.resize(n + 1);
  }
}

bool QEMHeap::greater(size_t i, size_t j) const {
  assert(!std::isnan(edges[keys[i]].error()));
  assert(!std::isnan(edges[keys[j]].error()));
  return edges[keys[i]].error() > edges[keys[j]].error();
}

void QEMHeap::exchange(size_t i, size_t j) {
  std::swap(handles[keys[i]], handles[keys[j]]);
  std::swap(keys[i], keys[j]);
}

bool QEMHeap::isMinHeap(size_t k) const {
  if (k > n) return true;
  size_t left = 2 * k;
  size_t right = 2 * k + 1;
  assert(!(left <= n && greater(k, left)));
  assert(!(right <= n && greater(k, right)));
  return isMinHeap(left) && isMinHeap(right);
}

void QEMHeap::swim(size_t k) {
  while (k > 1 && greater(k / 2, k)) {
    exchange(k, k / 2);
    k = k / 2;
  }
}

void QEMHeap::sink(size_t k) {
  assert(k >= 1);
  while (2 * k <= n) {
    size_t j = 2 * k;
    if (j < n && greater(j, j + 1)) ++j;
    if (!greater(k, j)) break;
    exchange(k, j);
    k = j;
  }
}

}  // namespace Internal
}  // namespace MeshSimpl
