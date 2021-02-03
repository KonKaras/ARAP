
#include <iostream>
#include <fstream>
#include <chrono>

#include "Eigen.h"
#include "SimpleMesh.h"
#include <Eigen/Sparse>
// #include "GUI.h"
#include "ArapDeformer.h"


using namespace std;

int main() {
	/***
	 *  SimpleMesh.m_vertices -> original positions
	 *  SimpleMesh.m_verticesPrime -> computed new position in deformed mesh
	 * 	ARAP.deform -> estimateVertices() löst das LGS Lx=b, wobei L=SimpleMesh.m_systemMatrix, x=pprime, b=b im code
	 * 
	 * 	test.off
	 * 	0--1--2--3--4  --> y
	 *  | /| /| /| /|
	 * 	5--6--7--8--9
	 * 	
	 * 	|
	 *  x
	 * 
	 ***/

	omp_set_num_threads(4);
	// Load the source and target mesh.
	const std::string filenameMesh = std::string("../data/bunny/cactus_small.off");
	
	vector<int> fixedPoints;

	fixedPoints.push_back(14);
	fixedPoints.push_back(13);
	fixedPoints.push_back(17);
	fixedPoints.push_back(259);
	fixedPoints.push_back(375);
	fixedPoints.push_back(388);
	fixedPoints.push_back(67);
	fixedPoints.push_back(586);

	int handleID = 586;
	Vector3f handleMoved(0, 0, 0.861985);

	SimpleMesh sourceMesh;
	if (!sourceMesh.loadMesh(filenameMesh)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
		return -1;
	}

	ArapDeformer deformer(&sourceMesh);
	
	auto t1 = std::chrono::high_resolution_clock::now();
	deformer.applyDeformation(fixedPoints, handleID, handleMoved, 3); 

	// fixedPoints.clear();
	// fixedPoints.push_back(0);
	// //fixedPoints.push_back(1);
	// // fixedPoints.push_back(2);
	// // fixedPoints.push_back(3);
	// fixedPoints.push_back(4);
	// // fixedPoints.push_back(5);
	// fixedPoints.push_back(7);
	// fixedPoints.push_back(10);
	// fixedPoints.push_back(14);
	// // fixedPoints.push_back(9);
	// handleID = 0;
	// handleMoved.x() = 0;
	// handleMoved.y() = 0;
	// handleMoved.z() = 2;

	// deformer.applyDeformation(fixedPoints, handleID, handleMoved, 3);


	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> eps = t2 - t1;
	std::cout << "Deformation completed in "<< eps.count() <<" seconds." << std::endl;

	// sourceMesh.writeMesh("../data/bunny/deformedMesh.off"); 

	return 0;
}
