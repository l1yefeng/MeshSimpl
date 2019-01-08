//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_QEM_HEAP_H
#define MESH_SIMPL_QEM_HEAP_H

#include "library.h"
#include <limits>
#include <boost/heap/binomial_heap.hpp>

namespace MeshSimpl
{

namespace Internal
{

typedef std::pair<idx, Edge*> QEMHeapNode;

struct QEMHeapCmp
{
    bool operator()(const QEMHeapNode& a, const QEMHeapNode& b) const
    {
        return a.second->error > b.second->error;
    }
};

class QEMHeap
{
public:
    // Construct a min-binary-heap with edge ecol errors as keys;
    // store a reference of the list of edges and store all handles
    explicit QEMHeap(std::vector<Edge>& edges) :edges(edges), min_heap(), handles(edges.size())
    {
        for (idx i = 0; i < edges.size(); ++i)
            if (edges[i].boundary_v != BOUNDARY_V::BOTH)
                handles[i] = min_heap.emplace(i, &edges[i]);
    }
    // Returns the edge id with minimum ecol error
    idx top() const { return min_heap.top().first; }
    // Remove the top edge from heap
    void pop() { min_heap.pop(); }
    // Fix the priority of an edge after the error value is modified
    void fix(idx e) { min_heap.update(handles[e]); }
    // Suppress this edge until it is, if ever, updated next time
    void penalize(idx e)
    {
        (*handles[e]).second->error = std::numeric_limits<double>::max();
        min_heap.decrease(handles[e]);
    }
    void erase(idx e) { min_heap.erase(handles[e]); }
    // Returns true if heap is empty
    bool empty() const { return min_heap.empty(); }
    // Returns the size of the heap, which should be the number of edges with boundary_v != BOTH
    size_t size() const { return min_heap.size(); }

private:
    const std::vector<Edge>& edges;
    boost::heap::binomial_heap<QEMHeapNode, boost::heap::compare<QEMHeapCmp>> min_heap;
    std::vector<boost::heap::binomial_heap<QEMHeapNode, boost::heap::compare<QEMHeapCmp>>::handle_type> handles;
};

}

}

#endif // MESH_SIMPL_QEM_HEAP_H
