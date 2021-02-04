#pragma once

#include <iostream>
#include <cstring>
#include <fstream>

#include "Eigen.h"
#include <unsupported/Eigen/src/MatrixFunctions/MatrixSquareRoot.h>
#include "SimpleMesh.h"
#include "ArapDeformer.h"

#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/unproject.h>

#include <list>


class GUI {
public:

	GUI(std::string filenameMesh, int iter, int get_weight_type, int get_estimation_type) {
		weight_type = get_weight_type;
		estimation_type = get_estimation_type;
		num_iterations = iter;
		displayMesh(filenameMesh);
	}

	std::set<int> staticFaces, staticFacesPreviousInit;
	std::set<int> handles, handlesPreviousInit;
	std::set<int> staticVertices, staticVerticesPreviousInit;
	int num_iterations;
	int weight_type;
	int estimation_type;

private:

	bool arapMode = false;
	bool handleSelectionMode = false;
	bool faceSelection = false;
	bool mouseDown = false;
	bool handleDown = false;
	bool vertexHit = false;
	bool arapInitialized = false;
	bool arapRunning = false;
	bool deformerInitiated = false;

	int currentMouseButton = 0;
	int currentMovingHandle = -1;
	int prevMovingHandle = -1;

	SimpleMesh sourceMesh;
	ArapDeformer deformer;
	std::string meshName;

	Eigen::MatrixXd vertices, colors, handleRep;
	Eigen::MatrixXi faces;

	//https://github.com/libigl/libigl/blob/master/tutorial/708_Picking/main.cpp

	void displayMesh(std::string filenameMesh) {
		//load mesh
		igl::readOFF(filenameMesh, vertices, faces);
		meshName = filenameMesh;

		//init white
		colors = Eigen::MatrixXd::Constant(faces.rows(), 3, 1);

		bool arapModeRef = arapMode;
		bool handleSelectionModeRef = handleSelectionMode;

		std::cout << "(Usage: [click]  Pick face on object)" << std::endl;
		std::cout << "(Usage: [press 1]  Select Static Vertices)" << std::endl;
		std::cout << "(Usage: [press 2]  Select Static Faces)" << std::endl;
		std::cout << "(Usage: [press 3]  Select Handles)" << std::endl;
		std::cout << "(Usage: [press 4]  Start ARAP Mode)" << std::endl;
		std::cout << "ARAP active" << std::endl;

		igl::opengl::glfw::Viewer viewer;

		bool paint = true;
		viewer.callback_mouse_down = [this, &paint](igl::opengl::glfw::Viewer& viewer, int button, int) -> bool
		{
			//checks if mouse key is pressed + we are not in arap mode + have selected a vertex with the mouse key press
			if (!arapMode) {
				currentMouseButton = button;
				mouseDown = true;

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
					if (handles.find(handleId) != handles.end()) {
						handleDown = true;
						currentMovingHandle = handleId;
						viewer.data().clear_points();
						set<int>::iterator itr;
						for (itr = staticVertices.begin(); itr != staticVertices.end(); itr++) {
							viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 0, 1));
						}
						RequestArapInit();
						return true;// DisplacementHandler(viewer);
					}
				}
			}
			return false;
		};

		viewer.callback_mouse_up = [this](igl::opengl::glfw::Viewer& viewer, int button, int) -> bool
		{
			if (button == currentMouseButton) {
				if (arapMode && handleDown) {
					set<int>::iterator itr;
					for (itr = handles.begin(); itr != handles.end(); itr++) {
						viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
					}
				}
				mouseDown = false;
				vertexHit = false;
				prevMovingHandle = currentMovingHandle;
				currentMovingHandle = -1;
				handleDown = false;
			}
			return false;
		};

		viewer.callback_mouse_move = [this](igl::opengl::glfw::Viewer& viewer, int, int) -> bool
		{
			//checks if mouse key is pressed + we are not in arap mode + have selected a vertex with the mouse key press
			if (mouseDown && vertexHit) {

				if (!arapMode) {
					return SelectionHandler(viewer, currentMouseButton);
				}
				else {
					if (handleDown) {
						return DisplacementHandler(viewer);
					}
				}
			}
			return false;
		};


		viewer.callback_key_down = [this](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) -> bool
		{
			//Fixed Vertex selection
			if (key == '1') {
				faceSelection = false;
				handleSelectionMode = false;
				arapMode = false;
				for each (int fid in staticFaces)
				{
					UpdateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
				}
				RepaintVertices(viewer);
				std::string out = "Static Vertices Selection active, ARAP Paused";
				std::cout << out << std::endl;
				return true;
			}
			if (key == '2') {
				faceSelection = true;
				handleSelectionMode = false;
				arapMode = false;
				for each (int fid in staticFaces)
				{
					UpdateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
				}
				RepaintVertices(viewer);
				std::string out = "Static Faces Selection active, ARAP Paused";
				std::cout << out << std::endl;
				return true;
			}
			//Toggle handle selection
			if (key == '3') {
				faceSelection = false;
				arapMode = false;
				for each (int fid in staticFaces)
				{
					UpdateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
				}
				RepaintVertices(viewer);
				handleSelectionMode = true;		
				std::cout << "Handle Selection Active, ARAP Paused" << std::endl;
			}
			//ARAP Mode
			if (key == '4') {
				faceSelection = false;
				handleSelectionMode = false;
				arapMode = true;
				for each (int fid in staticFaces)
				{
					UpdateColor(fid, Eigen::Vector3d(0, 0, 1), viewer);
				}
				//Repaint Vertices in blue
				viewer.data().clear_points();
				set<int>::iterator itr;
				for (itr = handles.begin(); itr != handles.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
				}
				for (itr = staticVertices.begin(); itr != staticVertices.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 0, 1));
				}
				std::string out = "ARAP active";
				std::cout << out << std::endl;
				return true;
			}
			return true;
		};

		//display mesh
		viewer.data().set_mesh(vertices, faces);

		//viewer.data().add_points(vertices, Eigen::RowVector3d(1, 0, 0));
		viewer.data().point_size = 10;
		viewer.data().set_colors(colors);
		viewer.data().show_lines = true;
		
		viewer.launch();
	}

	void PerformARAP(Eigen::Vector3f handlePos) {
		//only run algorithm if initialized and not already running
		if (arapInitialized && !arapRunning) {
			arapRunning = true;
			deformer.applyDeformation(GetStaticVerticesFromFaces(), currentMovingHandle, handlePos.cast<double>(), num_iterations); // Hier passiert die flipflop optimization mit 3 iterationen
			std::vector<Vertex> deformedVertices = deformer.m_mesh.getDeformedVertices();
			sourceMesh = deformer.m_mesh;
			cout << "GUI[handleID] is " << deformedVertices[currentMovingHandle].position.x() << "," << deformedVertices[currentMovingHandle].position.y() << "," << deformedVertices[currentMovingHandle].position.z() << endl;
			//MatrixXd deformedVerticesMat(deformedVertices.size(), 3);
			for (int i = 0; i < vertices.rows(); i++) {
				vertices.row(i) = deformedVertices[i].position.cast<double>();
			}
			arapRunning = false;
			//vertices = deformedVerticesMat;
			//std::cout << vertices.row(currentMovingHandle) << std::endl;
		}
	}

	void RequestArapInit() {
		//check if fixed or handle vertices have been added or removed -> re-init structs
		if (handles.size() != 0 && (staticFaces != staticFacesPreviousInit || staticVertices != staticVerticesPreviousInit || handles != handlesPreviousInit) || prevMovingHandle != currentMovingHandle) {
			arapInitialized = false;
			//prepare UI structs
			std::vector<int> staticsAsVector = GetStaticVerticesFromFaces();//(staticFaces.size());

			//std::cout << handlesAsVector.size() << std::endl;
			
			//std::cout << "handles done" << std::endl;
			//TODO adapt arap for multiple handles

			sourceMesh = SimpleMesh();
			
			if (!sourceMesh.loadMeshFromGUI(vertices, faces, staticsAsVector)) { // in loadMesh() finden wichtige vorberechnungen statt
				std::cout << "Mesh file wasn't read successfully at location: " << meshName << std::endl;
				return;
			}
			if (!deformerInitiated) {
				deformer = ArapDeformer(&sourceMesh, weight_type, estimation_type);
				deformerInitiated = true;
			}
			staticFacesPreviousInit = staticFaces;
			handlesPreviousInit = handles;
			staticVerticesPreviousInit = staticVertices;
			arapInitialized = true;
		}
	}

	std::vector<int> GetStaticVerticesFromFaces() {
		std::set<int> staticVerticesLocal;
		staticVerticesLocal.clear();
		staticVerticesLocal = staticVertices;
		//bool prevHandleIsNowStatic = false;
		for (int face : staticFaces) {
			for (int i = 0; i < 3; i++) {
				staticVerticesLocal.insert(faces.row(face)(i));
			}
		}
		//if the static faces do not contain the prev handle -> make prev handle non-static!
		//if (!prevHandleIsNowStatic && staticVertices.find(prevMovingHandle) != staticVertices.end()) staticVertices.erase(prevMovingHandle);

		if (staticVerticesLocal.find(currentMovingHandle) != staticVerticesLocal.end()) staticVerticesLocal.erase(currentMovingHandle);

		std::vector<int> staticVerticesAsVector(staticVerticesLocal.size());
		std::copy(staticVerticesLocal.begin(), staticVerticesLocal.end(), staticVerticesAsVector.begin());
		staticVerticesAsVector.push_back(currentMovingHandle);
		return staticVerticesAsVector;
	}

	bool DisplacementHandler(igl::opengl::glfw::Viewer& viewer) {
		double x = viewer.current_mouse_x;
		double y = viewer.core().viewport(3) - viewer.current_mouse_y;
		
		Eigen::Vector3f handlePos = vertices.row(currentMovingHandle).cast<float>();
		//Convert depth to view
		Eigen::Vector3d projection = igl::project(handlePos, viewer.core().view, viewer.core().proj, viewer.core().viewport).cast<double>();

		//Convert mouse position into world position
		Eigen::Vector3f worldPos = igl::unproject(Eigen::Vector3f(x, y, (float)projection.z()), viewer.core().view, viewer.core().proj, viewer.core().viewport);

		///vertices.row(currentMovingHandle) = worldPos.transpose();//+= diff;

		//TODO Send Data to ARAP
		PerformARAP(worldPos);
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
			/*auto& toSelect = handles;
			if (!handleSelectionMode) {
				if (faceSelection) {
					toSelect = staticFaces;
				}
				else {
					toSelect = staticVertices;
				}
			}*/
			//auto& toSelect = handleSelectionMode ? handles : staticFaces;
			if(!faceSelection) fid = GetClosestVertexIdFromBC(fid, bc);
			if (mouseID == 0) {
				if (handleSelectionMode) {
					Select(fid, bc, handles, viewer);
				}
				else {
					if (faceSelection) {
						Select(fid, bc, staticFaces, viewer);
					}
					else {
						Select(fid, bc, staticVertices, viewer);
					}
				}
			}
			else if(mouseID == 2){
				if (handleSelectionMode) {
					Unselect(fid, bc, handles, viewer);
				}
				else {
					if (faceSelection) {
						Unselect(fid, bc, staticFaces, viewer);
					}
					else {
						Unselect(fid, bc, staticVertices, viewer);
					}
				}
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
			if (!faceSelection) {
				viewer.data().clear_points();
				set<int>::iterator itr;
				for (itr = handles.begin(); itr != handles.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
				}
				for (itr = staticVertices.begin(); itr != staticVertices.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(1, 0, 0));
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

		if (!faceSelection) {
			if (handles.find(fid) == handles.end() && staticVertices.find(fid) == staticVertices.end()) {
				int selectedVertex = fid;//GetClosestVertexIdFromBC(fid, bc);
				toSelect.insert(selectedVertex);
				viewer.data().add_points(vertices.row(selectedVertex), newColor.transpose()); //Unsure about color vector
				//checks that vertex is hit
				vertexHit = true;
			}
		}
		else {
			//static faces are defined
			if (toSelect.find(fid) == toSelect.end()) {
				//select face
				toSelect.insert(fid);
				UpdateColor(fid, newColor, viewer);
				//std::cout << fid << std::endl;
			}
		}
	}

	void RepaintVertices(igl::opengl::glfw::Viewer& viewer) {
		viewer.data().clear_points();
		set<int>::iterator itr;
		for (itr = handles.begin(); itr != handles.end(); itr++) {
			viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
		}
		for (itr = staticVertices.begin(); itr != staticVertices.end(); itr++) {
			viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(1, 0, 0));
		}
	}

	
	int GetClosestVertexIdFromBC(int fid, Eigen::Vector3f& bc) {
		
		Eigen::Vector3i triangleVertices = faces.row(fid);
		//std::cout << triangleVertices << std::endl;
		
		Eigen::MatrixXd closestPoints (3,3);

		for (int i = 0; i < 3; i++) {
			closestPoints.row(i) = vertices.row(triangleVertices(i));
		}
		//std::cout << "closest vertices' coordinates" << std::endl;
		//std::cout << closestPoints << std::endl;
		
	
		Eigen::Vector3d sqdistances;
		Eigen::Vector3d queryPoint = bc(0) * closestPoints.row(0) + bc(1) * closestPoints.row(1)+ bc(2) * closestPoints.row(2);
		//std::cout << "Query" << std::endl;
		//std::cout << queryPoint << std::endl;

		
		//sqdistances = igl::point_mesh_squared_distance(Eigen::MatrixXd(queryPoint(0), queryPoint(1), queryPoint(2)), vertices, faces, sqdistances, closestPoints);
		for (int i = 0; i < 3; i++) {
			Eigen::Vector3d diff = closestPoints.row(i) - queryPoint.transpose();
			sqdistances(i) = std::sqrt(diff.dot(diff));
		}
		//std::cout << "Distances" << std::endl;
		//std::cout << sqdistances << std::endl;
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
		//std::cout << "smallest: " + std::to_string(min) + " for vertex " + std::to_string(minId);
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