#include <iostream>
#include <fstream>

#include "Eigen.h"
#include "VirtualSensor.h"
#include "SimpleMesh.h"
#include "ICPOptimizer.h"
#include "ARAPOptimizer.h"
#include "PointCloud.h"
using namespace std;


<<<<<<< HEAD

int main() {
	// Load the source and target mesh.
	const std::string filenameMesh = std::string("/home/nnrthmr/Code/MotionCaptring&3dScanning_WS20/FinalProject/arap/data/bunny/bunny.off");

	//Just for the moment some dummy points and handles
	std::vector<Vector3f> fixedPoints;
	fixedPoints.push_back(Vector3f(-0.0639191f, 0.179114f, -0.0588715f)); // right ear
	fixedPoints.push_back(Vector3f(0.0590575f, 0.066407f, 0.00686641f)); // tail
	fixedPoints.push_back(Vector3f(-0.0789843f, 0.13256f, 0.0519517f)); // mouth

	Vector3f handle(-0.0106867f, 0.179756f, -0.0283248f); // left ear

=======
int main() {
	// Load the source and target mesh.
	const std::string filenameMesh = std::string("/home/nnrthmr/Code/MotionCaptring&3dScanning_WS20/FinalProject/arap/data/bunny/bunny_part1.off");
>>>>>>> master

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

<<<<<<< HEAD
	// 1. Mesh load() : does all initializing and computes weight matrix and laplacian
	// 2. update laplacian if fixed points changed
	// 3. 

=======
>>>>>>> master
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
