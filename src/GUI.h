#pragma once

#include <iostream>
#include <cstring>
#include <fstream>

#include "Eigen.h"

#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

class GUI {
public:
	GUI() {

	}

	static void display(std::string filenameMesh) {
		//load mesh
		Eigen::MatrixXd vertices;
		Eigen::MatrixXi faces;

		igl::readOFF(filenameMesh, vertices, faces);

		//display mesh

		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(vertices, faces);
		viewer.launch();
	}
};