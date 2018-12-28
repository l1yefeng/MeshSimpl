#include <library.h>
#include <readOBJ.h>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " OBJFILE [STRENGTH]\n\n"
                  << "  OBJFILE   Path of input .obj file\n"
                  << "  STRENGTH  Simplification strength; The default is 0.5;\n"
                  << "            A strength of 0.8 will keep 20% vertices in the output\n"
                  << std::endl;
        return 0;
    }

    float strength = 0.5f;
    if (argc >= 3) {
        strength = std::stof(argv[2]);
        if (!(strength > 0 && strength < 1)) {
            std::cout << "Error: STRENGTH must be between 0 and 1\n"
                      << "  Use '" << argv[0] << "' without argument to view usage\n"
                      << std::endl;
            return 1;
        }
    }

    std::cout << "Received input: mesh data location: " << argv[1] << "\n"
              << "                simplification strength: " << strength << std::endl;

    // read obj file
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<unsigned int>> indices;
    igl::readOBJ(argv[1], vertices, indices);

    MeshSimpl::simplify(vertices, indices, strength);

    return 0;
}
