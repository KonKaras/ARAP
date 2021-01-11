#pragma once

#include <iostream>
#include <cstring>
#include <fstream>

#include "Eigen.h"

#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <list>

class GUI {
public:

	GUI(std::string filenameMesh) {
		display(filenameMesh);
	}

	std::set<int> selectedFaces;

	//https://github.com/libigl/libigl/blob/master/tutorial/708_Picking/main.cpp
	void display(std::string filenameMesh) {
		//load mesh

		Eigen::MatrixXd vertices, colors;
		Eigen::MatrixXi faces;

		auto& selectedFacesRef = selectedFaces;

		igl::readOFF(filenameMesh, vertices, faces);

		//init white
		colors = Eigen::MatrixXd::Constant(faces.rows(), 3, 1);

		bool arapMode = false;

		std::cout << "(Usage: [click]  Pick face on object)" << std::endl;
		std::cout << "(Usage: [press 1]  Toggle mode)" << std::endl;
		std::cout << "ARAP active" << std::endl;

		igl::opengl::glfw::Viewer viewer;

		viewer.callback_mouse_down = [this, &arapMode, &vertices, &faces, &colors, &selectedFacesRef] (igl::opengl::glfw::Viewer& viewer, int, int) -> bool
		{
			int fid;
			Eigen::Vector3f bc;
			// Cast a ray in the view direction starting from the mouse position
			double x = viewer.current_mouse_x;
			double y = viewer.core().viewport(3) - viewer.current_mouse_y;
			if (!arapMode) {
				if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
					viewer.core().proj, viewer.core().viewport, vertices, faces, fid, bc))
				{
					if (selectedFacesRef.find(fid) == selectedFacesRef.end()) {
						//select face
						selectedFacesRef.insert(fid);
						UpdateColor(fid, Eigen::Vector3d(1,0,0), colors, viewer);
						return true;
					}
					else {
						//unselect face
						selectedFacesRef.erase(fid);
						UpdateColor(fid, Eigen::Vector3d(1, 1, 1), colors, viewer);
						return true;
					}
				}
			}
			return false;
		};
		viewer.callback_key_down = [this, &arapMode, &colors, &selectedFacesRef](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) -> bool
		{
			if (key == '1') {
				arapMode = !arapMode;
				if (arapMode) {
					for each (int fid in selectedFacesRef)
					{
						UpdateColor(fid, Eigen::Vector3d(0, 0, 1), colors, viewer);
					}
				}
				else {
					for each (int fid in selectedFacesRef)
					{
						UpdateColor(fid, Eigen::Vector3d(1, 0, 0), colors, viewer);
					}
				}

				std::string out = arapMode ? "ARAP active" : "Selection active";
				std::cout << out << std::endl;
				return true;
			}
			return true;
		};

		//display mesh
		viewer.data().set_mesh(vertices, faces);
		viewer.data().set_colors(colors);
		viewer.data().show_lines = true;
		viewer.launch();
	}


private:
	void UpdateColor(int faceID, Eigen::Vector3d newColor, Eigen::MatrixXd& colors, igl::opengl::glfw::Viewer& viewer) {
		colors.row(faceID) << newColor(0), newColor(1), newColor(2);
		viewer.data().set_colors(colors);
	}
};