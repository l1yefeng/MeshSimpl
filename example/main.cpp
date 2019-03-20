#include "clipp.h"
#include "simplify.hpp"
#include <chrono>
#include <cstdlib>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace clipp;

int main(int argc, char* argv[]) {
    string in, out;
    bool fix_boundary;
    string weighting;
    float strength;

    auto cli = (value("input", in).doc("input .obj file"),
                option("-o", "--output") & value("path", out).doc("output file path"),
                option("-f", "--fix-boundary")
                    .set(fix_boundary)
                    .doc("do not move vertices on boundary"),
                option("-w", "--weighting") &
                    word("strategy", weighting).doc("one of [area, uniform]"),
                option("-s", "--strength") &
                    number("ratio", strength).doc("0.8 means remove 80% vertices"));

    if (!parse(argc, argv, cli)) {
        cout << make_man_page(cli, argv[0]);
        return 1;
    }

    MeshSimpl::SimplifyOptions options;
    options.strength = strength;
    options.fix_boundary = fix_boundary;
    if (weighting == "area")
        options.weighting = MeshSimpl::BY_AREA;

    // read obj file
    vector<vector<double>> vertices, out_vertices;
    vector<vector<unsigned int>> indices, out_indices;
    igl::readOBJ(in, vertices, indices);

    try {
        const auto before = chrono::steady_clock::now();
        auto res = MeshSimpl::simplify(vertices, indices, options);
        const auto after = chrono::steady_clock::now();
        const long duration =
            chrono::duration_cast<chrono::milliseconds>(after - before).count();
        std::cout << "[INFO] Simplification completed (" << duration << " milliseconds)"
                  << std::endl;
        out_vertices = res.first;
        out_indices = res.second;
    } catch (char const* exception) {
        cerr << exception << endl;
        return 1;
    }

    Eigen::MatrixXd V(out_vertices.size(), 3);
    for (int i = 0; i < V.rows(); ++i)
        V.row(i) << out_vertices[i][0], out_vertices[i][1], out_vertices[i][2];
    Eigen::MatrixXi F(out_indices.size(), 3);
    for (int i = 0; i < F.rows(); ++i)
        F.row(i) << out_indices[i][0], out_indices[i][1], out_indices[i][2];

    if (out.empty()) {
        // view
        Eigen::MatrixXd OV(vertices.size(), 3);
        for (int i = 0; i < OV.rows(); ++i)
            OV.row(i) << vertices[i][0], vertices[i][1], vertices[i][2];
        Eigen::MatrixXi OF(indices.size(), 3);
        for (int i = 0; i < OF.rows(); ++i)
            OF.row(i) << indices[i][0], indices[i][1], indices[i][2];

        bool show_simplified = true;
        const auto& key_down_cb = [&](igl::opengl::glfw::Viewer& viewer,
                                      unsigned char key, int mod) -> bool {
            switch (key) {
            case 's':
            case 'S':
                show_simplified = !show_simplified;
                break;
            default:
                return false;
            }
            viewer.data().clear();
            if (show_simplified)
                viewer.data().set_mesh(V, F);
            else
                viewer.data().set_mesh(OV, OF);
            return true;
        };

        igl::opengl::glfw::Viewer viewer;
        viewer.callback_key_down = key_down_cb;
        viewer.data().set_mesh(V, F);
        viewer.launch();
    } else {
        // write
        igl::writeOBJ(out, V, F);
    }
}
