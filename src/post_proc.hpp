//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_POST_PROC_HPP
#define LIB_MESH_SIMPL_POST_PROC_HPP

#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

// Remove the vertices and indices that are marked deleted, and reduce the vector size
void compact_data(const std::vector<bool>& deleted_vertex,
                  const std::vector<bool>& deleted_face, V& vertices, F& indices);

} // namespace Internal
} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_POST_PROC_HPP
