//
// Created by nickl on 1/3/19.
//

#ifndef MESH_SIMPL_WRITE_OBJ_H
#define MESH_SIMPL_WRITE_OBJ_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

void write_obj(const std::string& filename, const std::vector<std::vector<double>>& vertices,
               const std::vector<std::vector<unsigned int>>& indices) {
    const int precision = 9;

    int coord_col_len = precision + 4;
    double coord = 0.0;
    for (const auto& v : vertices)
        coord = std::max({coord, std::abs(v[0]), std::abs(v[1]), std::abs(v[2])});
    for (unsigned sz = std::floor(coord); sz > 0; sz /= 10, ++coord_col_len)
        ;

    int v_col_len = 1;
    for (unsigned sz = vertices.size(); sz > 0; sz /= 10, ++v_col_len)
        ;

    std::ofstream ofs(filename);

    ofs << "#" << std::endl
        << "# number of vertices: " << vertices.size() << std::endl
        << "# number of faces:    " << indices.size() << std::endl
        << "#" << std::endl
        << std::endl;

    ofs << std::fixed << std::setprecision(precision);
    for (const auto& v : vertices) {
        ofs << "v";
        for (const double x : v)
            ofs << std::right << std::setw(coord_col_len) << x;
        ofs << std::endl;
    }
    ofs << std::endl;
    for (const auto& f : indices) {
        ofs << "f";
        for (const unsigned int v : f)
            ofs << std::right << std::setw(v_col_len) << v + 1;
        ofs << std::endl;
    }
    ofs << std::endl;

    ofs.close();
}

#endif // MESH_SIMPL_WRITE_OBJ_H
