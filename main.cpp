#include <iostream>
#include <library.h>
#include <readOBJ.h>

int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " OBJFILE [STRENGTH]\n\n"
			  << "  OBJFILE   Path of input .obj file\n"
			  << "  STRENGTH  Simplification strength; The default is 0.5;\n"
			  << "            A strength of 0.8 will keep 20% vertices in the output\n"
			  << std::endl;
		return 1;
	}

	float strength = 0.5f;
	if (argc >= 3) {
		strength = std::stof(argv[2]);
	}

	std::cout << "Received input: mesh data location: " << argv[1] << "\n"
		  << "                simplification strength: " << strength << std::endl;

	// read obj file
	Eigen::MatrixXd vertices;
	Eigen::MatrixXi indices;
	igl::readOBJ(argv[1], vertices, indices);
	// print some stats of mesh
	std::cout << "Mesh info: #V = " << vertices.rows() << "\n"
		  << "           #F = " << indices.rows() << std::endl;

	MeshSimpl::hello();

	return 0;
}
