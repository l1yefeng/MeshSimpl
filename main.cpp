#include <library.h>
#include <readOBJ.h>

using namespace std;

void write_obj(const string& filename,
               const vector<vector<double>>& vertices,
               const vector<vector<unsigned int>>& indices);

int main(int argc, char* argv[])
{
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
    }
    catch (const std::invalid_argument& e) {
        cerr << "Error: STRENGTH cannot be recognized" << endl;
        return 1;
    }
    if (!(strength > 0 && strength < 1)) {
        cout << "Error: STRENGTH must be between 0 and 1" << endl
             << "  Use '" << argv[0] << "' without argument to view usage" << endl
             << endl;
        return 1;
    }

    cout << "Received input: mesh data location: " << argv[1] << endl
         << "                simplification strength: " << strength << endl;

    // read obj file
    vector<vector<double>> vertices;
    vector<vector<unsigned int>> indices;
    igl::readOBJ(argv[1], vertices, indices);

    const auto res = MeshSimpl::simplify(vertices, indices, strength);
    write_obj(argv[2], res.first, res.second);

    return 0;
}

void write_obj(const string& filename,
               const vector<vector<double>>& vertices,
               const vector<vector<unsigned int>>& indices)
{
    ofstream ofs(filename);

    ofs << "#" << endl
        << "# author: nickl" << endl
        << "#" << endl
        << "# number of vertices: " << vertices.size() << endl
        << "# number of faces: " << indices.size() << endl
        << "#" << endl
        << endl;

    for (const auto& v : vertices) {
        ofs << "v";
        for (const double x : v)
            ofs << " " << x;
        ofs << endl;
    }

    ofs << endl;

    for (const auto& f : indices) {
        ofs << "f";
        for (const unsigned int v : f)
            ofs << " " << v;
        ofs << endl;
    }

    ofs.close();
}

