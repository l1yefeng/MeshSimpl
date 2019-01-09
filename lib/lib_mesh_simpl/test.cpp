//
// Created by nickl on 12/29/18.
//

#define CATCH_CONFIG_MAIN

#include "util.h"
#include "ecol.h"
#include "qem_heap.h"
#include "post_proc.h"
#include <catch2/catch.hpp>

using namespace std;
using namespace MeshSimpl::Internal;

TEST_CASE("Quadric error should be zero when v is on the plane", "[q_error]")
{
    Quadric q = make_quadric({0, 1, 0}, -0.5);
    double err = q_error(q, std::vector<double>{1, 0.5, 1});
    REQUIRE(err == 0.0);
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
    compact_data(v_del, f_del, V, F);

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

TEST_CASE("QEM heap should behave normally", "[QEMHeap]")
{
    vector<Edge> edges(10);
    for (int i = 0; i < 10; ++i) {
        edges[i].error = 1.1*((i+5)%10);
        edges[i].boundary_v = BOUNDARY_V::NONE;
    }
    Edge edge_should_ignore{};
    edge_should_ignore.boundary_v = BOUNDARY_V::BOTH;
    edges.push_back(edge_should_ignore);

    QEMHeap heap(edges);
    REQUIRE(heap.size() == 10);

    SECTION("constructor") {
        array<unsigned int, 10> e_results{5, 6, 7, 8, 9, 0, 1, 2, 3, 4};
        for (auto e : e_results) {
            REQUIRE(heap.top() == e);
            heap.pop();
        }
        REQUIRE(heap.empty());
    }

    SECTION("size()") {
        REQUIRE(heap.size() == 10);
    }

    SECTION("fix()") {
        edges[0].error = -1.0;
        heap.fix(&edges[0], false);
        edges[7].error = 8.0;
        heap.fix(&edges[7], true);
        array<unsigned int, 10> e_results{0, 5, 6, 8, 9, 1, 2, 7, 3, 4};
        for (auto e : e_results) {
            REQUIRE(heap.top() == e);
            heap.pop();
        }
    }

    SECTION("penalize()") {
        heap.penalize(0);
        heap.penalize(7);
        array<unsigned int, 8> e_results{5, 6, 8, 9, 1, 2, 3, 4};
        for (auto e : e_results) {
            REQUIRE(heap.top() == e);
            heap.pop();
        }
    }

    SECTION("erase()") {
        heap.erase(2);
        heap.erase(0);
        heap.erase(7);
        REQUIRE(heap.size() == 7);
        array<unsigned int, 7> e_results{5, 6, 8, 9, 1, 3, 4};
        for (auto e : e_results) {
            REQUIRE(heap.top() == e);
            heap.pop();
        }
    }
}
