//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_QEM_HEAP_H
#define MESH_SIMPL_QEM_HEAP_H

#include "types.h"
#include <limits>
#include <vector>

namespace MeshSimpl
{

namespace Internal
{

// Reference: https://algs4.cs.princeton.edu/24pq/MinPQ.java
class QEMHeap
{
public:
    // Construct a min-binary-heap with edge ecol errors as keys;
    // store a reference of the list of edges and store all handles
    explicit QEMHeap(std::vector<Edge>& edges) :edges(edges), keys(edges.size()+1),
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
    // Returns the edge id with minimum ecol error
    idx top() const { return keys[1]; };
    // Remove the top edge from heap
    void pop()
    {
        exchange(1, n--);
        sink(1);
        keys.resize(n+1);
        assert(is_min_heap());
    }
    // Fix the priority of an edge after the error value is modified;
    // Param sinking should be true if the node will have lower priority after fix-up
    void fix(Edge* const ptr, bool sinking)
    {
        idx e = static_cast<idx>(ptr-edges.data());
        size_t k = handles[e];
        if (sinking)
            sink(k);
        else
            swim(k);
    }
    // Suppress this edge until it is, if ever, updated next time
    void penalize(idx e)
    {
        edges[e].error = std::numeric_limits<double>::max();
        sink(handles[e]);
    }
    void erase(idx e) {
        const double error_prev = edges[e].error;
        size_t k = handles[e];
        exchange(k, n--);
        if (edges[keys[k]].error > error_prev)
            sink(k);
        else
            swim(k);
    }
    // Returns true if heap is empty
    bool empty() const { return n == 0; };
    // Returns the size of the heap, which should be the number of edges with boundary_v != BOTH
    size_t size() const { return n; };

private:
    std::vector<idx> keys;
    std::vector<Edge>& edges;
    std::vector<size_t> handles;
    size_t n;

    // Compare function: larger error --> lower priority
    bool greater(size_t i, size_t j) const { return edges[keys[i]].error > edges[keys[j]].error; }
    // Helper function: the sole function that modifies handles and keys data
    void exchange(size_t i, size_t j)
    {
        std::swap(handles[keys[i]], handles[keys[j]]);
        std::swap(keys[i], keys[j]);
    }
    // For assertion purposes
    bool is_min_heap() const { return is_min_heap(1); }
    bool is_min_heap(size_t k) const
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
    // Adjust the priority of node k: direction is towards higher priority
    void swim(size_t k)
    {
        while (k > 1 && greater(k/2, k)) {
            exchange(k, k/2);
            k = k/2;
        }
    }
    // Adjust the priority of node k: direction is towards lower priority
    void sink(size_t k)
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
};

}

}

#endif // MESH_SIMPL_QEM_HEAP_H
