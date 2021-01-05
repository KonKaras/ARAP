#include <iostream>
#include <fstream>

#include "Eigen.h"
#include "VirtualSensor.h"
#include "SimpleMesh.h"
#include "ICPOptimizer.h"
#include "ARAPOptimizer.h"
#include "PointCloud.h"
using namespace std;


int main() {
	// Load the source and target mesh.
	const std::string filenameMesh = std::string("/home/nnrthmr/Code/MotionCaptring&3dScanning_WS20/FinalProject/arap/data/bunny/bunny_part1.off");

	SimpleMesh sourceMesh;
	if (!sourceMesh.loadMesh(filenameMesh)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
		return -1;
	}

	// ARAPOptimizer *optimizer = new ARAPOptimizer();
	// optimizer->setNumIterations(10);

	// PointCloud mesh{ sourceMesh };
	// Matrix4f estimatedPose = optimizer->deform();

	std::cout << "ARAP done." << std::endl;

	// delete optimizer;

	/*
	d = Deformer(filename)
	d.read_file()
	d.build_weight_matrix()
	if len(selection_filename) > 0:
		d.read_deformation_file(deformation_file)
		d.read_selection_file(selection_filename)
	d.calculate_laplacian_matrix()
	d.precompute_p_i()
	print("Precomputation time ", time.time() - t)
	t = time.time()
	d.apply_deformation(iterations)
	print("Total iteration time", time.time() - t)
	d.output_s_prime_to_file()
	d.show_graph()
	*/

	return 0;
}
