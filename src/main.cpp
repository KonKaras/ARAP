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
	const std::string filenameMesh = std::string("../../data/bunny/test.off");

	vector<int> fixedPoints;
	// fixedPoints.push_back(6); 
	// fixedPoints.push_back(7); 
	// fixedPoints.push_back(8);  

	// int handleID = 140; // left ear
	// Vector4f handleMoved(2*0.037177, 2*0.110644, 2*0.018811, 1); // left ear moved to new position (2*original position)

	fixedPoints.push_back(2);
	fixedPoints.push_back(7);
	int handleID = 0;
	Vector4f handleMoved(0,0, 2, 1);


	SimpleMesh sourceMesh;
	if (!sourceMesh.loadMesh(filenameMesh, fixedPoints, handleID)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
		return -1;
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	applyDeformation(sourceMesh, handleID, handleMoved, 1);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> eps = t2 - t1;
	std::cout << "Deformation completed in "<< eps.count() <<" seconds." << std::endl;

	sourceMesh.writeMesh("../../data/bunny/deformedMesh.off");

	return 0;
}
