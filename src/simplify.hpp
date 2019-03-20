//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_SIMPLIFY_HPP
#define LIB_MESH_SIMPL_SIMPLIFY_HPP

#include "types.hpp"

namespace MeshSimpl {

struct SimplifyOptions {
    // simplifies until vertex count is 1-strength of the original,
    // only accept value in range [0, 1)
    float strength = 0.5f;

    // weight the quadrics by triangle area, i.e., Q becomes Q * weight
    // a larger weight makes the computed error larger thus "later" to modify
    // during the iterations of edge collapse operations.
    //  - UNIFORM: no weighting
    //  - BY_AREA: larger face -> larger error
    WEIGHTING weighting = UNIFORM;

    // when "fix_boundary" is true, we completely do not collapse anything on boundary;
    // otherwise, we add a "constraint plane" that is perpendicular to boundary face
    // to increase the quadric/error of boundary vertices
    // ref: Simplifying Surfaces with Color and Texture using Quadric Error Metrics
    bool fix_boundary = false;
};

void options_validation(const SimplifyOptions& options);

// Mesh simplification main method.
// TODO: Ignoring the 4th and the following values (if exist) in vertices.
std::pair<V, F> simplify(const V& vertices, const F& indices,
                         const SimplifyOptions& options = {});

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_SIMPLIFY_HPP
