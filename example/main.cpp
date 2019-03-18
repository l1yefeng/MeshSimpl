#include <chrono>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <simplify.hpp>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " OBJFILE OUTPUT [STRENGTH]" << endl
             << endl
             << "  OBJFILE   Path of input .obj file" << endl
             << "  OUTPUT    Path of output data" << endl
             << "  STRENGTH  Simplification strength; The default is 0.5;" << endl
             << "            A strength of 0.8 will keep 20% vertices in the output"
             << endl
             << endl;
        return argc == 1 ? 0 : 1;
    }

    float strength;
    try {
        strength = stof(argv[3]);
    } catch (const std::invalid_argument& e) {
        cerr << "Error: STRENGTH cannot be recognized" << endl;
        return 1;
    }
    if (!(strength > 0 && strength < 1)) {
        cout << "Error: STRENGTH must be between 0 and 1" << endl
             << "  Use '" << argv[0] << "' without argument to view usage" << endl
             << endl;
        return 1;
    }

    // read obj file
    vector<vector<double>> vertices;
    vector<vector<unsigned int>> indices;
    igl::readOBJ(argv[1], vertices, indices);

    MeshSimpl::SimplifyOptions options;
    options.strength = strength;
    options.weighting = MeshSimpl::BY_AREA;
    options.fix_boundary = false;

    const auto before = chrono::steady_clock::now();

    try {
        auto res = MeshSimpl::simplify(vertices, indices, options);
        const auto after = chrono::steady_clock::now();
        const long duration =
            chrono::duration_cast<chrono::milliseconds>(after - before).count();
        std::cout << "[INFO] Simplification completed (" << duration << " milliseconds)"
                  << std::endl;

        Eigen::MatrixXd V(res.first.size(), 3);
        for (int i = 0; i < V.rows(); ++i)
            V.row(i) << res.first[i][0], res.first[i][1], res.first[i][2];
        Eigen::MatrixXi F(res.second.size(), 3);
        for (int i = 0; i < F.rows(); ++i)
            F.row(i) << res.second[i][0], res.second[i][1], res.second[i][2];

        igl::writeOBJ(argv[2], V, F);
    } catch (char const* exception) {
        cerr << exception << endl;
        return 1;
    }
}