#include "measure.h"
#include "write_obj.h"
#include <readOBJ.h>
#include <simplify.h>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " OBJFILE OUTPUT [STRENGTH]" << endl
             << endl
             << "  OBJFILE   Path of input .obj file" << endl
             << "  OUTPUT    Path of output data" << endl
             << "  STRENGTH  Simplification strength; The default is 0.5;" << endl
             << "            A strength of 0.8 will keep 20% vertices in the output" << endl
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
    options.debug = true;
//    Measure measure;
    try {
        const auto res = MeshSimpl::simplify(vertices, indices, options);
//        long duration = measure.stop();
//        std::cout << "[INFO] Simplification completed (" << duration << " milliseconds)"
//                  << std::endl;
        write_obj(argv[2], res.first, res.second);
    } catch (char const* exception) {
        cerr << exception << endl;
        return 1;
    }
}
