#include <Windows.h>
#include <iostream>
#include <fstream>
#include <chrono>

#include "Eigen.h"
#include "VirtualSensor.h"
#include "SimpleMesh.h"
#include "ICPOptimizer.h"
#include "Arap.h"
#include "PointCloud.h"
#include "GUI.h"


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

	// Load the source and target mesh.
	const std::string filenameMesh = std::string("../data/bunny/test.off");
	bool debug = false;
	if (!debug) {
		GUI* gui = new GUI(filenameMesh, 3);
	}
	else {
		
		vector<int> fixedPoints; //hier sammeln wir die vertices die sich nciht bewegen durfen
		fixedPoints.push_back(0);
		//fixedPoints.push_back(1);
		fixedPoints.push_back(2);
		fixedPoints.push_back(3);
		fixedPoints.push_back(4);
		fixedPoints.push_back(5);
		fixedPoints.push_back(6);
		fixedPoints.push_back(7);
		fixedPoints.push_back(8);
		fixedPoints.push_back(9);
		int handleID = 2;
		Vector3f handleMoved(-3, 2, 0);


		SimpleMesh sourceMesh;
		if (!sourceMesh.loadMesh(filenameMesh, fixedPoints)) { // in loadMesh() finden wichtige vorberechnungen statt
			std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
			return -1;
		}

		auto t1 = std::chrono::high_resolution_clock::now();
		applyDeformation(&sourceMesh, handleID, handleMoved, 3); // Hier passiert die flipflop optimization mit 3 iterationen
		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> eps = t2 - t1;
		std::cout << "Deformation completed in "<< eps.count() <<" seconds." << std::endl;

		sourceMesh.writeMesh("../data/bunny/deformedMesh.off");
		
	}
	return 0;
}
