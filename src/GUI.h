#pragma once

#include <iostream>
#include <cstring>
#include <fstream>

#include "Eigen.h"
#include <unsupported/Eigen/src/MatrixFunctions/MatrixSquareRoot.h>

#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/point_mesh_squared_distance.h>

#include <list>


class GUI {
public:

	GUI(std::string filenameMesh) {
		displayMesh(filenameMesh);
	}

	std::set<int> staticFaces;
	std::set<int> faceHandles;

private:

	bool arapMode = false;
	bool handleSelectionMode = false;

	Eigen::MatrixXd vertices, colors;
	Eigen::MatrixXi faces;

	//https://github.com/libigl/libigl/blob/master/tutorial/708_Picking/main.cpp

	void displayMesh(std::string filenameMesh) {
		//load mesh

		igl::readOFF(filenameMesh, vertices, faces);

		//init white
		colors = Eigen::MatrixXd::Constant(faces.rows(), 3, 1);

		bool arapModeRef = arapMode;
		bool handleSelectionModeRef = handleSelectionMode;

		std::cout << "(Usage: [click]  Pick face on object)" << std::endl;
		std::cout << "(Usage: [press 1]  Toggle between ARAP and Selection)" << std::endl;
		std::cout << "(Usage: [press 2]  Select Handles)" << std::endl;
		std::cout << "ARAP active" << std::endl;

		igl::opengl::glfw::Viewer viewer;

		viewer.callback_mouse_down = [this](igl::opengl::glfw::Viewer& viewer, int, int) -> bool
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
					ClickerHandler(fid, bc, viewer);
					return true;
				}
				return false;
			}
			//this for free rotation during selection
			return false;
		};


		viewer.callback_key_down = [this](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) -> bool
		{
			//Toggle ARAP
			if (key == '1') {
				handleSelectionMode = false;
				arapMode = !arapMode;
				if (arapMode) {
					for each (int fid in staticFaces)
					{
						UpdateColor(fid, Eigen::Vector3d(0, 0, 1), viewer);
					}
				}
				else {
					for each (int fid in staticFaces)
					{
						UpdateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
					}
				}

				std::string out = arapMode ? "ARAP active" : "Static Face Selection active";
				std::cout << out << std::endl;
				return true;
			}
			//Toggle handle selection
			if (key == '2') {
				arapMode = false;
				for each (int fid in staticFaces)
				{
					UpdateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
				}
				handleSelectionMode = true;		
				std::cout << "Handle Selection Active, ARAP Paused" << std::endl;
			}
			return true;
		};

		//display mesh
		viewer.data().set_mesh(vertices, faces);

		//viewer.data().add_points(vertices, Eigen::RowVector3d(1, 0, 0));

		viewer.data().set_colors(colors);
		viewer.data().show_lines = true;
		viewer.data().show_face_labels = true;
		
		viewer.launch();
	}

	bool ClickerHandler(int fid, Eigen::Vector3f& bc, igl::opengl::glfw::Viewer& viewer)
	{
		auto& faces = handleSelectionMode ? faceHandles : staticFaces;
		auto newColor = handleSelectionMode ? Vector3d(0, 1, 0) : Vector3d(1, 0, 0);

		if (handleSelectionMode) {
			int selectedVertex = GetClosestVertexIdFromBC(fid, bc);
			std::cout << selectedVertex << std::endl;
			viewer.data().add_points(vertices.row(selectedVertex), Eigen::RowVector3d(0,1,0));
		}
		else {
			//static faces are defined
			if (faces.find(fid) == faces.end()) {
				//select face
				faces.insert(fid);
				UpdateColor(fid, newColor, viewer);
				//std::cout << fid << std::endl;
			}
			else {
				//unselect face
				faces.erase(fid);
				UpdateColor(fid, Eigen::Vector3d(1, 1, 1), viewer);
			}

			//remove from other set if necessary, a handle cannot be static and a static face cannot be a handle
			auto& otherFaces = handleSelectionMode ? staticFaces : faceHandles;
			if (otherFaces.find(fid) != otherFaces.end()) {
				otherFaces.erase(fid);
				//std::cout << "erased face " + std::to_string(fid) + " from " + (handleSelectionMode ? "staticFaces" : "faceHandles") << std::endl;
			}
		}
		return true;
	}

	
	int GetClosestVertexIdFromBC(int fid, Eigen::Vector3f& bc) {
		
		Eigen::Vector3i triangleVertices = faces.row(fid);
		std::cout << triangleVertices << std::endl;
		
		Eigen::MatrixXd closestPoints (3,3);

		for (int i = 0; i < 3; i++) {
			closestPoints.row(i) = vertices.row(triangleVertices(i));
		}
		std::cout << "closest vertices' coordinates" << std::endl;
		std::cout << closestPoints << std::endl;
		
	
		Eigen::Vector3d sqdistances;
		Eigen::Vector3d queryPoint = bc(0) * closestPoints.row(0) + bc(1) * closestPoints.row(0)+ bc(2) * closestPoints.row(2);
		std::cout << "Query" << std::endl;
		std::cout << queryPoint << std::endl;

		
		//sqdistances = igl::point_mesh_squared_distance(Eigen::MatrixXd(queryPoint(0), queryPoint(1), queryPoint(2)), vertices, faces, sqdistances, closestPoints);
		for (int i = 0; i < 3; i++) {
			Eigen::Vector3d diff = closestPoints.row(i) - queryPoint.transpose();
			sqdistances(i) = std::sqrt(diff.dot(diff));
		}
		std::cout << "Distances" << std::endl;
		std::cout << sqdistances << std::endl;
		float min = 100;
		int minId = 0;
		float current;
		for (int i = 0; i < 3; i++) {
			current = sqdistances(i);
			if (min > current) {
				min = current;
				minId = i;
			}
		}
		std::cout << "smallest: " + std::to_string(min) + " for vertex " + std::to_string(minId);
		return triangleVertices(minId);
		
	}

	void UpdateColor(int faceID, Eigen::Vector3d newColor, igl::opengl::glfw::Viewer& viewer) {
		colors.row(faceID) << newColor(0), newColor(1), newColor(2);
		viewer.data().set_colors(colors);
	}
};


//TODO make a nice menu with dropdown or buttons for modes, possibly brush size etc
class Menu {
public:
	Menu(igl::opengl::glfw::Viewer viewer)
	{
		this->viewer = viewer;
	}
private:
	igl::opengl::glfw::Viewer viewer;
	
};