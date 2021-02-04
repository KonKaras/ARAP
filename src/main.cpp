#define NOMINMAX
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <chrono>

#include "Eigen.h"
#include "SimpleMesh.h"
#include <Eigen/Sparse>
#include "GUI.h"
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

	// Load the source and target mesh.
	bool isInvalid = true;
	std::string filenameMesh, filename;
	do {
		isInvalid = true;
		filenameMesh = std::string("../data/bunny/");
		std::cout << "Insert filename: " << std::endl;
		std::cin >> filename;
		filenameMesh.append(filename);
		filenameMesh.append(".off");

		std::ifstream file(filenameMesh);
		if (file.is_open()) {
			isInvalid = false;
		}
		if (isInvalid) {
			std::cout << "Mesh file does not exist." << std::endl;
			std::cin.clear();
			std::cin.ignore();
		}
	} while (isInvalid);

	int weight_type, estimation_type;
	do {
		isInvalid = true;
		std::cout << "Choose type of weighting function (0 = uniform, 1 = constant, 2 = cotangent)" << std::endl;
		std::cin >> weight_type;

		if (!std::cin.fail()) {
			if (weight_type >= 0 && weight_type <= 2) {
				isInvalid = false;
			}
		}
		if (isInvalid) {
			std::cout << "Invalid input!" << std::endl;
			std::cin.clear();
			std::cin.ignore();
		}
	} while (isInvalid);

	do {
		isInvalid = true;
		std::cout << "Choose type of matrix decomposition (0 = sparse QR, 1 = sparse LU, 2 = non sparce matrices)" << std::endl;
		std::cin >> estimation_type;

		if (!std::cin.fail()) {
			if (estimation_type >= 0 && estimation_type <= 2) {
				isInvalid = false;
			}
		}
		if (isInvalid) {
			std::cout << "Invalid input!" << std::endl;
			std::cin.clear();
			std::cin.ignore();
		}
	} while (isInvalid);

	bool debug = false;
	if (!debug) {
		GUI* gui = new GUI(filenameMesh, 100, weight_type, estimation_type);
    }
    else{
		vector<int> fixedPoints;
		for (int i = 0; i <= 20; i++) {
			fixedPoints.push_back(i);
			fixedPoints.push_back(420+i);
		}
		fixedPoints.push_back(199);

		SimpleMesh sourceMesh;
		if (!sourceMesh.loadMesh(filenameMesh)) {
			std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
			return -1;
		}

		ArapDeformer deformer(&sourceMesh, weight_type, estimation_type);

		int handleID = 199;
		Vector3d handleMoved = sourceMesh.getVertex(handleID) + Vector3d(0,0,0.5);

		auto t1 = std::chrono::high_resolution_clock::now();
		//deformer.initDeformation(fixedPoints);
		deformer.applyDeformation(fixedPoints, handleID, handleMoved, 100);
		deformer.m_mesh.writeMesh("../data/bunny/deformedMesh1.off");
		/*
		fixedPoints.clear();
		fixedPoints.push_back(0);
		//fixedPoints.push_back(1);
		fixedPoints.push_back(2);
		fixedPoints.push_back(3);
		fixedPoints.push_back(4);
		// fixedPoints.push_back(5);
		fixedPoints.push_back(7);
		//fixedPoints.push_back(10);
		//fixedPoints.push_back(14);
		// fixedPoints.push_back(9);
		handleID = 7;
		handleMoved = Vector3d(1, 2, 2);
		//deformer.initDeformation(fixedPoints);
		deformer.applyDeformation(fixedPoints, handleID, handleMoved, 3);
		*/
		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> eps = t2 - t1;
		std::cout << "Deformation completed in "<< eps.count() <<" seconds." << std::endl;
	}
	return 0;
}
