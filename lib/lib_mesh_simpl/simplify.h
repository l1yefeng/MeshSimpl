//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_SIMPLIFY_H
#define LIB_MESH_SIMPL_SIMPLIFY_H

#include "types.h"

namespace MeshSimpl {

struct SimplifyOptions {
    float strength = 0.5f;      // simplifies until vertex count is 1-strength of the original
    bool weight_by_face = true; // weight quadrics by triangle area
    bool debug = false;         // print some debugging info to standard output
};

void options_validation(const SimplifyOptions& options);

// Mesh simplification main method. Simplify given mesh until remaining number of vertices/faces
// is (1-strength) of the original. Returns output vertices and indices as in inputs.
// TODO: Ignoring the 4th and the following values (if exist) in vertices.
std::pair<V, F> simplify(const V& vertices, const F& indices, const SimplifyOptions& options = {});

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_SIMPLIFY_H
