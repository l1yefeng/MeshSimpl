//
// Created by nickl on 1/8/19.
//

#include "qem_heap.h"
#include "ecol.h"
#include <cassert>

namespace MeshSimpl
{

namespace Internal
{

QEMHeap::QEMHeap(std::vector<Edge>& edges) :keys(edges.size()+1), edges(edges),
                                            handles(edges.size()), n(0)
{
    for (idx i = 0; i < edges.size(); ++i)
        if (edges[i].boundary_v != BOUNDARY_V::BOTH)
            keys[handles[i] = ++n] = i;

    keys.resize(n);
    for (size_t k = n/2; k >= 1; --k)
        sink(k);
    assert(is_min_heap());
}

void QEMHeap::pop()
{
    exchange(1, n--);
    sink(1);
    keys.resize(n+1);
    assert(is_min_heap());
}

void QEMHeap::fix(Edge* const ptr, bool sinking)
{
    auto e = static_cast<idx>(ptr-edges.data());
    size_t k = handles[e];
    if (sinking)
        sink(k);
    else
        swim(k);
}

void QEMHeap::penalize(idx e)
{
    edges[e].error = std::numeric_limits<double>::max();
    sink(handles[e]);
}

void QEMHeap::erase(idx e)
{
    const double error_prev = edges[e].error;
    size_t k = handles[e];
    exchange(k, n--);
    if (edges[keys[k]].error > error_prev)
        sink(k);
    else
        swim(k);
}

void QEMHeap::reset_edge_errors(const V& vertices, const std::vector<Quadric>& quadrics)
{
    assert(keys.size() == n+1);
    // re-compute error and collapse center for living edges
    for (size_t k = 1; k < keys.size(); ++k) {
        auto& edge = edges[keys[k]];
        set_edge_error(vertices, quadrics, edge);
    }

    // fix priority queue
    for (size_t k = n/2; k >= 1; --k)
        sink(k);
    assert(is_min_heap());
}

void QEMHeap::exchange(size_t i, size_t j)
{
    std::swap(handles[keys[i]], handles[keys[j]]);
    std::swap(keys[i], keys[j]);
}

bool QEMHeap::is_min_heap(size_t k) const
{
    if (k > n)
        return true;
    size_t left = 2*k;
    size_t right = 2*k+1;
    if (left <= n && greater(k, left))
        return false;
    if (right <= n && greater(k, right))
        return false;
    return is_min_heap(left) && is_min_heap(right);
}

void QEMHeap::swim(size_t k)
{
    while (k > 1 && greater(k/2, k)) {
        exchange(k, k/2);
        k = k/2;
    }
}

void QEMHeap::sink(size_t k)
{
    while (2*k <= n) {
        size_t j = 2*k;
        if (j < n && greater(j, j+1))
            ++j;
        if (!greater(k, j))
            break;
        exchange(k, j);
        k = j;
    }
}

}

}
