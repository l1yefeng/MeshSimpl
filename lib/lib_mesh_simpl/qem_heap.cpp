//
// Created by nickl on 1/8/19.
//

#include "qem_heap.h"
#include "ecol.h"
#include <cassert>
#include <limits>

namespace MeshSimpl {

namespace Internal {

QEMHeap::QEMHeap(E& edges) : keys(edges.size() + 1), edges(edges), handles(edges.size()), n(0) {
    for (idx i = 0; i < edges.size(); ++i)
        if (edges[i].boundary_v != BOUNDARY_V::BOTH)
            keys[handles[i] = ++n] = i;

    keys.resize(n + 1);
    prioritize_all();
}

void QEMHeap::pop() {
    exchange(1, n--);
    sink(1);
    keys.resize(n + 1);
}

void QEMHeap::fix(Edge* const ptr, bool sinking) {
    auto e = static_cast<idx>(ptr - edges.data());
    size_t k = handles[e];
    if (sinking)
        sink(k);
    else
        swim(k);
}

void QEMHeap::penalize(idx e) {
    edges[e].error = std::numeric_limits<double>::max();
    sink(handles[e]);
}

void QEMHeap::erase(idx e) {
    const double error_prev = edges[e].error;
    size_t k = handles[e];
    exchange(k, n--);
    if (edges[keys[k]].error > error_prev)
        sink(k);
    else
        swim(k);
    keys.resize(n + 1);
}

void QEMHeap::exchange(size_t i, size_t j) {
    std::swap(handles[keys[i]], handles[keys[j]]);
    std::swap(keys[i], keys[j]);
}

bool QEMHeap::is_min_heap(size_t k) const {
    if (k > n)
        return true;
    size_t left = 2 * k;
    size_t right = 2 * k + 1;
    if (left <= n && greater(k, left))
        return false;
    if (right <= n && greater(k, right))
        return false;
    return is_min_heap(left) && is_min_heap(right);
}

void QEMHeap::swim(size_t k) {
    while (k > 1 && greater(k / 2, k)) {
        exchange(k, k / 2);
        k = k / 2;
    }
}

void QEMHeap::sink(size_t k) {
    while (2 * k <= n) {
        size_t j = 2 * k;
        if (j < n && greater(j, j + 1))
            ++j;
        if (!greater(k, j))
            break;
        exchange(k, j);
        k = j;
    }
}
void QEMHeap::prioritize_all() {
    for (size_t k = n / 2; k >= 1; --k)
        sink(k);
    assert(is_min_heap());
}

} // namespace Internal

} // namespace MeshSimpl
