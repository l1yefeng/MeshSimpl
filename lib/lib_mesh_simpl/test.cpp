//
// Created by nickl on 12/29/18.
//

#define CATCH_CONFIG_MAIN

#include "util.h"
#include "library.h"
#include <algorithm>
#include <catch2/catch.hpp>

using namespace std;
using namespace MeshSimpl::Internal;

TEST_CASE("Quadric error should be zero when v is on the plane", "[q_error]")
{
    Quadric q = make_quadric({0, 1, 0}, -0.5);
    double err = q_error(q, {1, 0.5, 1}, false);
    REQUIRE(err == 0.0);
}

TEST_CASE("The edge with smallest error should be popped", "[build_min_heap]")
{
    vector<Edge> edges(3);
    for (auto& edge : edges)
        edge.boundary_v = BOUNDARY_V::NONE;
    edges[0].error = 0.1;
    edges[1].error = 2.2;
    edges[2].error = 0.7;
    edges.push_back({});
    edges[3].boundary_v = BOUNDARY_V::BOTH;
    auto heap = build_min_heap(edges);

    REQUIRE(heap.size() == 3);
    REQUIRE(heap.top() == 0);
    heap.pop();
    REQUIRE(heap.top() == 2);
}

TEST_CASE("Deleted values should be gone", "[compact_data]")
{
    vector<vector<double>> V{
        {-0.500000, -0.500000, 0.500000},
        {0.500000,  -0.500000, 0.500000},
        {-0.500000, 0.500000,  0.500000},
        {0.500000,  0.500000,  0.500000},
        {-0.500000, 0.500000,  -0.500000},
        {0.500000,  0.500000,  -0.500000},
        {-0.500000, -0.500000, -0.500000},
        {0.500000,  -0.500000, -0.500000},
    };
    vector<vector<unsigned int>> F{
        {0, 1, 2},
        {2, 1, 3},
        {2, 3, 4},
        {4, 3, 5},
        {4, 5, 6},
        {6, 5, 7},
        {6, 7, 0},
        {0, 7, 1},
        {1, 7, 3},
        {3, 7, 5},
        {6, 0, 4},
        {4, 0, 2},
    };
    vector<bool> v_del(V.size(), false);
    v_del[0] = true;
    vector<bool> f_del(F.size(), false);
    for (int i : {0, 6, 7, 10, 11})
        f_del[i] = true;
    compact_data(V, F, v_del, f_del);

    SECTION("Vertices") {
        REQUIRE(V.size() == 7);
    }

    SECTION("Indices") {
        REQUIRE(F.size() == 7);
        auto it = find_if(F.begin(), F.end(), [&](const vector<unsigned int>& face)->bool {
            return !(face[0] < V.size() && face[1] < V.size() && face[2] < V.size());
        });
        REQUIRE(it == F.end());
    }
}
