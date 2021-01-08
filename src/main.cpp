#include <iostream>
#include <fstream>
#include <chrono>

#include "Eigen.h"
#include "VirtualSensor.h"
#include "SimpleMesh.h"
#include "ICPOptimizer.h"
#include "Arap.h"
#include "PointCloud.h"
using namespace std;


int main() {
	// Load the source and target mesh.
	const std::string filenameMesh = std::string("../../data/bunny/bunny.off");

	vector<int> fixedPoints;
	fixedPoints.push_back(6); 
	fixedPoints.push_back(7); 
	fixedPoints.push_back(8);  

	int handle = 140; // left ear
	Vector4f handleMoved(2*0.037177, 2*0.110644, 2*0.018811, 1); // left ear moved

	SimpleMesh sourceMesh;
	if (!sourceMesh.loadMesh(filenameMesh, fixedPoints, handle)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
		return -1;
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	applyDeformation(sourceMesh, 1);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> eps = t2 - t1;
	std::cout << "Deformation completed in "<< eps.count() <<" seconds." << std::endl;

	//sourceMesh.writeMesh("../../data/bunny/deformedMesh.off");


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
	*/

	return 0;
}
