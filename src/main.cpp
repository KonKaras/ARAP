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
	const std::string filenameMesh = std::string("/home/nnrthmr/Code/MotionCaptring&3dScanning_WS20/FinalProject/arap/data/bunny/bunny.off");

	vector<Vector4f> fixedPoints;
	fixedPoints.push_back(Vector4f(0.0060412, 0.1249436, 0.03265289, 1)); 
	fixedPoints.push_back(Vector4f(-0.01346903, 0.1629936, -0.01200002, 1)); 
	fixedPoints.push_back(Vector4f(-0.03439324, 0.1723669, -0.0009821299, 1));  

	Vector4f handle(0.037177, 0.110644, 0.018811, 1); // left ear

	SimpleMesh sourceMesh;
	if (!sourceMesh.loadMesh(filenameMesh, fixedPoints)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
		return -1;
	}

	// ARAPOptimizer *optimizer = new ARAPOptimizer();
	// optimizer->setNumIterations(10);

	// PointCloud mesh{ sourceMesh };
	// Matrix4f estimatedPose = optimizer->deform();

	std::cout << "ARAP done." << std::endl;

	// delete optimizer;


	// 1. Mesh load() : does all initializing and computes weight matrix and laplacian
	// 2. update laplacian if fixed points changed
	// 3. 

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
