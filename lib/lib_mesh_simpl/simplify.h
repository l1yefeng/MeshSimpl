//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_SIMPLIFY_H
#define LIB_MESH_SIMPL_SIMPLIFY_H

#include "types.h"

namespace MeshSimpl {

struct SimplifyOptions {
    // simplifies until vertex count is 1-strength of the original
    float strength = 0;
    // weight quadrics by triangle area
    WEIGHTING weighting = BY_AREA;
    // force boundary to be fixed
    bool fix_boundary = false;
};

void options_validation(const SimplifyOptions& options);

// Mesh simplification main method.
// TODO: Ignoring the 4th and the following values (if exist) in vertices.
std::pair<V, F> simplify(const V& vertices, const F& indices,
                         const SimplifyOptions& options = {});

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_SIMPLIFY_H
