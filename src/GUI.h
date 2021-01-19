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
#include <igl/unproject.h>

#include <list>


class GUI {
public:

	GUI(std::string filenameMesh) {
		displayMesh(filenameMesh);
	}

	std::set<int> staticFaces;
	std::set<int> handles;

private:

	bool arapMode = false;
	bool handleSelectionMode = false;
	bool mouseDown = false;
	bool vertexHit = false;

	int currentMouseButton = 0;
	int currentMovingHandle = 0;

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

		bool paint = true;
		viewer.callback_mouse_down = [this, &paint](igl::opengl::glfw::Viewer& viewer, int button, int) -> bool
		{
			//checks if mouse key is pressed + we are not in arap mode + have selected a vertex with the mouse key press
			if (!arapMode) {
				currentMouseButton = button;
				mouseDown = true;
				/*
				int fid;
				Eigen::Vector3f bc;
				// Cast a ray in the view direction starting from the mouse position
				double x = viewer.current_mouse_x;
				double y = viewer.core().viewport(3) - viewer.current_mouse_y;
				if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
					viewer.core().proj, viewer.core().viewport, vertices, faces, fid, bc))
				{
					vertexHit = true;
					ClickerHandler(fid, bc, viewer);
					return true;
				}
				return false;
				*/
				return SelectionHandler(viewer, button);
			}
			else {
				currentMouseButton = button;
				mouseDown = true;

				int fid;
				Eigen::Vector3f bc;

				double x = viewer.current_mouse_x;
				double y = viewer.core().viewport(3) - viewer.current_mouse_y;

				if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
					viewer.core().proj, viewer.core().viewport, vertices, faces, fid, bc))
				{
					vertexHit = true;
					int handleId = GetClosestVertexIdFromBC(fid, bc);
					if (handles.find(handleId) == handles.end()) {
						currentMovingHandle = handleId;
						return DisplacementHandler(viewer);
					}
				}
			}
			return false;
		};

		viewer.callback_mouse_up = [this](igl::opengl::glfw::Viewer& viewer, int button, int) -> bool
		{
			if (button == currentMouseButton) {
				mouseDown = false;
				vertexHit = false;
			}
			return false;
		};

		viewer.callback_mouse_move = [this](igl::opengl::glfw::Viewer& viewer, int, int) -> bool
		{
			//checks if mouse key is pressed + we are not in arap mode + have selected a vertex with the mouse key press
			if (mouseDown && vertexHit) {
				/*
				int fid;
				Eigen::Vector3f bc;
				// Cast a ray in the view direction starting from the mouse position
				double x = viewer.current_mouse_x;
				double y = viewer.core().viewport(3) - viewer.current_mouse_y;
				if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
					viewer.core().proj, viewer.core().viewport, vertices, faces, fid, bc))
				{
					ClickerHandler(fid, bc, viewer);
					return true;
				}
				*/
				if (!arapMode) {
					return SelectionHandler(viewer, currentMouseButton);
				}
				else {
					return DisplacementHandler(viewer);
				}
			}
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

	bool DisplacementHandler(igl::opengl::glfw::Viewer& viewer) {
		int fid;
		Eigen::Vector3f bc;

		double x = viewer.current_mouse_x;
		double y = viewer.core().viewport(3) - viewer.current_mouse_y;

		//Transform mouse pos into world pos
		Eigen::Vector2d mousePos = 1 / vertices.row(currentMovingHandle).z() * Eigen::Vector2d(x,y);
		//Eigen::Matrix4d projection = viewer.
		Eigen::Vector3d worldPos = igl::unproject(Eigen::Vector3d(x, y, 1 / vertices.row(currentMovingHandle).z()), viewer.core().view, viewer.core().proj, viewer.core().viewport);

		Eigen::Vector3d diff = vertices.row(currentMovingHandle) - Eigen::Vector3d(x, y, vertices.row(currentMovingHandle).z()).transpose();
		vertices.row(currentMovingHandle) += diff;

		//TODO Send Data to ARAP
		//ARAP does stuff
		//repaint
		viewer.data().set_mesh(vertices, faces);
		return true;
	}

	bool SelectionHandler(igl::opengl::glfw::Viewer& viewer, int& mouseID) {
		int fid;
		Eigen::Vector3f bc;
		// Cast a ray in the view direction starting from the mouse position
		double x = viewer.current_mouse_x;
		double y = viewer.core().viewport(3) - viewer.current_mouse_y;
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
			viewer.core().proj, viewer.core().viewport, vertices, faces, fid, bc))
		{
			vertexHit = true;
			auto& toSelect = handleSelectionMode ? handles : staticFaces;
			if(handleSelectionMode) fid = GetClosestVertexIdFromBC(fid, bc);
			if (mouseID == 0) {
				Select(fid, bc, toSelect, viewer);
			}
			else if(mouseID == 2){
				Unselect(fid, bc, toSelect, viewer);
			}
			else {
				return false;
			}
			return true;
		}
		return false;
	}

	void Unselect(int fid, Eigen::Vector3f& bc, std::set<int>& toSelect, igl::opengl::glfw::Viewer& viewer) {
		if (toSelect.find(fid) != toSelect.end()) {
			toSelect.erase(fid);
			//repaint
			if (handleSelectionMode) {
				viewer.data().clear_points();
				set<int>::iterator itr;
				for (itr = handles.begin(); itr != handles.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
				}
			}
			else {
				UpdateColor(fid, Eigen::Vector3d(1, 1, 1), viewer);
			}
		}
	}

	void Select(int fid, Eigen::Vector3f& bc, std::set<int>& toSelect, igl::opengl::glfw::Viewer& viewer)
	{
		auto newColor = handleSelectionMode ? Vector3d(0, 1, 0) : Vector3d(1, 0, 0);

		if (handleSelectionMode) {
			int selectedVertex = fid;//GetClosestVertexIdFromBC(fid, bc);
			toSelect.insert(selectedVertex);
			viewer.data().add_points(vertices.row(selectedVertex), Eigen::RowVector3d(0,1,0));
			//checks that vertex is hit
			vertexHit = true;
		}
		else {
			//static faces are defined
			if (toSelect.find(fid) == toSelect.end()) {
				//select face
				toSelect.insert(fid);
				UpdateColor(fid, newColor, viewer);
				//std::cout << fid << std::endl;
			}
			/*
			else {
				//unselect face
				toSelect.erase(fid);
				UpdateColor(fid, Eigen::Vector3d(1, 1, 1), viewer);
			}
			
			//remove from other set if necessary, a handle cannot be static and a static face cannot be a handle
			auto& otherFaces = handleSelectionMode ? staticFaces : handles;
			if (otherFaces.find(fid) != otherFaces.end()) {
				otherFaces.erase(fid);
				//std::cout << "erased face " + std::to_string(fid) + " from " + (handleSelectionMode ? "staticFaces" : "faceHandles") << std::endl;
			}
			*/
		}
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
		Eigen::Vector3d queryPoint = bc(0) * closestPoints.row(0) + bc(1) * closestPoints.row(1)+ bc(2) * closestPoints.row(2);
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