//
// Created by nickl on 12/29/18.
//

#define CATCH_CONFIG_MAIN

#include "util.h"
#include <algorithm>
#include <catch2/catch.hpp>

using namespace std;

TEST_CASE("Quadric can be sumed", "[Quadric arithmetic]")
{
    typedef MeshSimpl::Internal::Quadric Q;

    Q q0{};
    Q q1{2, 3, 2, 1, 2, 0.1, 2.3, 2.3, 2.3, 7};
    Q q2{3, 4, 3, 2, 3, 1.1, 3.3, 3.3, 3.3, 4};

    SECTION("operator+=") {
        q0 += q1;
        REQUIRE(equal(q0.begin(), q0.end(), q1.begin()));
    }

    SECTION("operator+") {
        Q res{};
        for (unsigned int i = 0; i < 10; ++i)
            res[i] = q1[i]+q2[i];
        q0 = q1+q2;
        REQUIRE(equal(q0.begin(), q0.end(), res.begin()));
    }
}

