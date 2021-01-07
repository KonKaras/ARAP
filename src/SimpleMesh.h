#pragma once

#include <iostream>
#include <fstream>

#include "Eigen.h"
using namespace std;

struct Vertex {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	// Position stored as 4 floats (4th component is supposed to be 1.0)
	Vector4f position;
	// Color stored as 4 unsigned char
	Vector4uc color;
};

struct Triangle {
	unsigned int idx0;
	unsigned int idx1;
	unsigned int idx2;

	Triangle() : idx0{ 0 }, idx1{ 0 }, idx2{ 0 } {}

	Triangle(unsigned int _idx0, unsigned int _idx1, unsigned int _idx2) :
		idx0(_idx0), idx1(_idx1), idx2(_idx2) {}

};


class SimpleMesh {
public:
	SimpleMesh() {}

	/**
	 * Constructs a mesh from the current color and depth image.
	 */
	SimpleMesh(VirtualSensor& sensor, const Matrix4f& cameraPose, float edgeThreshold = 0.01f) {
		// Get ptr to the current depth frame.
		// Depth is stored in row major (get dimensions via sensor.GetDepthImageWidth() / GetDepthImageHeight()).
		float* depthMap = sensor.getDepth();
		// Get ptr to the current color frame.
		// Color is stored as RGBX in row major (4 byte values per pixel, get dimensions via sensor.GetColorImageWidth() / GetColorImageHeight()).
		BYTE* colorMap = sensor.getColorRGBX();

		// Get depth intrinsics.
		Matrix3f depthIntrinsics = sensor.getDepthIntrinsics();
		float fovX = depthIntrinsics(0, 0);
		float fovY = depthIntrinsics(1, 1);
		float cX = depthIntrinsics(0, 2);
		float cY = depthIntrinsics(1, 2);

		// Compute inverse depth extrinsics.
		Matrix4f depthExtrinsicsInv = sensor.getDepthExtrinsics().inverse();

		// Compute inverse camera pose (mapping from camera CS to world CS).
		Matrix4f cameraPoseInverse = cameraPose.inverse();

		// Compute vertices with back-projection.
		m_vertices.resize(sensor.getDepthImageWidth() * sensor.getDepthImageHeight());
		// For every pixel row.
		for (unsigned int v = 0; v < sensor.getDepthImageHeight(); ++v) {
			// For every pixel in a row.
			for (unsigned int u = 0; u < sensor.getDepthImageWidth(); ++u) {
				unsigned int idx = v*sensor.getDepthImageWidth() + u; // linearized index
				float depth = depthMap[idx];
				if (depth == MINF) {
					m_vertices[idx].position = Vector4f(MINF, MINF, MINF, MINF);
					m_vertices[idx].color = Vector4uc(0, 0, 0, 0);
				}
				else {
					// Back-projection and tranformation to world space.
					m_vertices[idx].position = cameraPoseInverse * depthExtrinsicsInv * Vector4f((u - cX) / fovX * depth, (v - cY) / fovY * depth, depth, 1.0f);

					// Project position to color map.
					Vector3f proj = sensor.getColorIntrinsics() * (sensor.getColorExtrinsics() * cameraPose * m_vertices[idx].position).block<3, 1>(0, 0);
					proj /= proj.z(); // dehomogenization
					unsigned int uCol = (unsigned int)std::floor(proj.x());
					unsigned int vCol = (unsigned int)std::floor(proj.y());
					if (uCol >= sensor.getColorImageWidth()) uCol = sensor.getColorImageWidth() - 1;
					if (vCol >= sensor.getColorImageHeight()) vCol = sensor.getColorImageHeight() - 1;
					unsigned int idxCol = vCol*sensor.getColorImageWidth() + uCol; // linearized index color
																					//unsigned int idxCol = idx; // linearized index color

					// Write color to vertex.
					m_vertices[idx].color = Vector4uc(colorMap[4 * idxCol + 0], colorMap[4 * idxCol + 1], colorMap[4 * idxCol + 2], colorMap[4 * idxCol + 3]);
				}
			}
		}

		// Compute triangles (faces).
		m_triangles.reserve((sensor.getDepthImageHeight() - 1) * (sensor.getDepthImageWidth() - 1) * 2);
		for (unsigned int i = 0; i < sensor.getDepthImageHeight() - 1; i++) {
			for (unsigned int j = 0; j < sensor.getDepthImageWidth() - 1; j++) {
				unsigned int i0 = i*sensor.getDepthImageWidth() + j;
				unsigned int i1 = (i + 1)*sensor.getDepthImageWidth() + j;
				unsigned int i2 = i*sensor.getDepthImageWidth() + j + 1;
				unsigned int i3 = (i + 1)*sensor.getDepthImageWidth() + j + 1;

				bool valid0 = m_vertices[i0].position.allFinite();
				bool valid1 = m_vertices[i1].position.allFinite();
				bool valid2 = m_vertices[i2].position.allFinite();
				bool valid3 = m_vertices[i3].position.allFinite();

				if (valid0 && valid1 && valid2) {
					float d0 = (m_vertices[i0].position - m_vertices[i1].position).norm();
					float d1 = (m_vertices[i0].position - m_vertices[i2].position).norm();
					float d2 = (m_vertices[i1].position - m_vertices[i2].position).norm();
					if (edgeThreshold > d0 && edgeThreshold > d1 && edgeThreshold > d2)
						addFace(i0, i1, i2);
				}
				if (valid1 && valid2 && valid3) {
					float d0 = (m_vertices[i3].position - m_vertices[i1].position).norm();
					float d1 = (m_vertices[i3].position - m_vertices[i2].position).norm();
					float d2 = (m_vertices[i1].position - m_vertices[i2].position).norm();
					if (edgeThreshold > d0 && edgeThreshold > d1 && edgeThreshold > d2)
						addFace(i1, i3, i2);
				}
			}
		}
	}

	void clear() {
		m_vertices.clear();
		m_triangles.clear();
	}

	unsigned int addVertex(Vertex& vertex) {
		unsigned int vId = (unsigned int)m_vertices.size();
		m_vertices.push_back(vertex);
		return vId;
	}

	unsigned int addFace(unsigned int idx0, unsigned int idx1, unsigned int idx2) {
		unsigned int fId = (unsigned int)m_triangles.size();
		Triangle triangle(idx0, idx1, idx2);
		m_triangles.push_back(triangle);
		return fId;
	}

	std::vector<Vertex>& getVertices() {
		return m_vertices;
	}

	const std::vector<Vertex>& getVertices() const {
		return m_vertices;
	}

	std::vector<Triangle>& getTriangles() {
		return m_triangles;
	}

	const std::vector<Triangle>& getTriangles() const {
		return m_triangles;
	}

	void transform(const Matrix4f& transformation) {
		for (Vertex& v : m_vertices) {
			v.position = transformation * v.position;
		}
	}

	bool loadMesh(const std::string& filename) {
		// Read off file (Important: Only .off files are supported).
		m_vertices.clear();
		m_verticesPrime.clear();
		m_triangles.clear();

		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cout << "Mesh file wasn't read successfully." << std::endl;
			return false;
		}

		// First line should say 'COFF'.
		char string1[5];
		file >> string1;

		// Read header.
		unsigned int numV = 0; //vertices
		unsigned int numP = 0; //faces
		unsigned int numE = 0; //edges
		file >> numV >> numP >> numE;

		m_vertices.reserve(numV);
		m_verticesPrime.reserve(numV);
		m_triangles.reserve(numP);

		m_numV = numV;

		// Read vertices.
		if (std::string(string1).compare("COFF") == 0) {
			// We have color information.
			for (unsigned int i = 0; i < numV; i++) {
				Vertex v;
				file >> v.position.x() >> v.position.y() >> v.position.z();
				v.position.w() = 1.f;
				// Colors are stored as integers. We need to convert them.
				Vector4i colorInt;
				file >> colorInt.x() >> colorInt.y() >> colorInt.z() >> colorInt.w();
				v.color = Vector4uc((unsigned char)colorInt.x(), (unsigned char)colorInt.y(), (unsigned char)colorInt.z(), (unsigned char)colorInt.w());
				m_vertices.push_back(v);
				m_verticesPrime.push_back(v);
			}
		}
		else if (std::string(string1).compare("OFF") == 0) {
			// We only have vertex information.
			for (unsigned int i = 0; i < numV; i++) {
				Vertex v;
				file >> v.position.x() >> v.position.y() >> v.position.z();
				v.position.w() = 1.f;
				v.color.x() = 0;
				v.color.y() = 0;
				v.color.z() = 0;
				v.color.w() = 255;
				m_vertices.push_back(v);
				m_verticesPrime.push_back(v);
			}
		}
		else {
			std::cout << "Incorrect mesh file type." << std::endl;
			return false;
		}


		//Create adjacency matrix
		m_neighborMatrix = MatrixXf::Zero(m_numV, m_numV);
		m_verticesToFaces = std::vector<std::vector<unsigned int>>(m_numV); //vtf[i] contains faces that contain vertex ID i
        m_cellRotations = std::vector<MatrixXf>(m_numV); // list of Rs per fan
		m_edgeMatrix = MatrixXf::Zero(m_numV,m_numV);

		// Read faces (i.e. triangles).
		for (unsigned int i = 0; i < numP; i++) {
			unsigned int num_vs;
			file >> num_vs;
			ASSERT(num_vs == 3 && "We can only read triangular mesh.");
			
			Triangle t;
			file >> t.idx0 >> t.idx1 >> t.idx2;
			m_triangles.push_back(t);

			//fill adjacency matrix
			addFaceToAdjacencyMatrix(t.idx0, t.idx1, t.idx2);  // 1 for neighbors, 0 else
			(m_verticesToFaces[t.idx0]).push_back(i); // list of facenumbers for each vertex
			(m_verticesToFaces[t.idx1]).push_back(i);
			(m_verticesToFaces[t.idx2]).push_back(i);
		}


		for(int i=0; i<numV;i++){
		 	m_edgeMatrix(i,i) = (m_neighborMatrix.rowwise().sum())[i];
			m_cellRotations[i] = MatrixXf::Zero(3,3);
		}

		cout << "numVertices: "<<numV<<endl;
		cout << "numFaces: "<<numP<<endl;
		cout << "numEdges: "<<numE<<endl;

		cout << m_edgeMatrix.rows() << " - " << m_edgeMatrix.cols()<< endl;
		cout << m_neighborMatrix.rows() << " - " << m_neighborMatrix.cols()<< endl;

		// m_laplaceMatrix = m_edgeMatrix - m_neighborMatrix; // TODO necessary?

		//buildWeightMatrix();
		// computeDistances();


		return true;
	}

	void addFaceToAdjacencyMatrix(int v1, int v2 ,int v3){
        m_neighborMatrix(v1, v2) = 1;
        m_neighborMatrix(v2, v1) = 1;
        m_neighborMatrix(v1, v3) = 1;
        m_neighborMatrix(v3, v1) = 1;
        m_neighborMatrix(v2, v3) = 1;
        m_neighborMatrix(v3, v2) = 1;
	}

	bool partOfFace(int i, int j, Triangle f){
		return (f.idx0 ==i && (f.idx1==j || f.idx2==j)) ||
				(f.idx1 ==i && (f.idx0==j || f.idx2==j)) ||
				(f.idx2 ==i && (f.idx1==j || f.idx0==j)) ;
	}

	vector<int> getNeighborsOf(int vertId){
        vector<int> neighbors;
        for(int i=0; i<m_numV; i++){
			if(m_neighborMatrix(vertId, i) == 1)
				neighbors.push_back(i);
		}
        return neighbors;
	}

	int getThirdFacePoint(int i, int j, Triangle f){
		vector<int> points;
		points.push_back(f.idx0);
		points.push_back(f.idx1);
		points.push_back(f.idx2);
		for(int idx : points){
			if(idx!=i && idx!=j){
				return idx;
			}
		}
		return -1;
	}

	void buildWeightMatrix(){

		cout<<"Generating Weight Matrix"<<endl;
        m_weightMatrix = MatrixXf(m_numV, m_numV);
		m_weightSum = MatrixXf(m_numV, m_numV);

        for (int vertex_id= 0; vertex_id < m_numV;vertex_id++){
			vector<int> neighbors = getNeighborsOf(vertex_id);
            for (int neighbor_id : neighbors)
                assignWeightForPair(vertex_id, neighbor_id);
		}

	}

	void assignWeightForPair(int i, int j){
		float weightIJ;
        if(m_weightMatrix(j, i) == 0) //If the opposite weight has not been computed, do so
            weightIJ = computeWeightForPair(i, j);
        else

            weightIJ = m_weightMatrix(j, i);

        m_weightSum(i, i) += (weightIJ * 0.5);
        m_weightSum(j, j) += (weightIJ * 0.5);
        m_weightMatrix(i, j) = weightIJ;
	}
	
	double angleBetweenVectors(Vector4f a, Vector4f b) {
		Vector3f a3(a.x(), a.y(), a.z());
		Vector3f b3(b.x(), b.y(), b.z());
		return std::atan2(a3.cross(b3).norm(), a3.dot(b3));
	}

    float computeWeightForPair(int i, int j){
        vector<Triangle> localFaces;
        for(unsigned int id : m_verticesToFaces[i]){
			Triangle f = m_triangles[id];
            if(partOfFace(i,j,f)) //If the face contains both I and J, add it
                localFaces.push_back(f);
		}
            

        // Either a normal face or a boundry edge, otherwise bad mesh
        assert(localFaces.size() <= 2);

        Vertex vertex_i = m_vertices[i];
        Vertex vertex_j = m_vertices[j];

        // Using the following weights: 0.5 * (cot(alpha) + cot(beta))

        float cot_theta_sum = 0.0f;
        for (Triangle f : localFaces){
			int other_vertex_id = getThirdFacePoint(i, j, f);
            Vertex vertex_o = m_vertices[other_vertex_id];

            float theta = angleBetweenVectors(vertex_i.position - vertex_o.position, vertex_j.position - vertex_o.position);
            cot_theta_sum += (1 / tan(theta));
		}
            
        return cot_theta_sum * 0.5;
	}

	void calculateLaplaceMatrix(){ // TODO so far this makes additional constraints for each fixed vertex, do we need this if we use ceres?
        
        m_laplaceMatrix = m_weightSum - m_weightMatrix;
        int num_fixedVertices = m_fixedVertices.size();

        int n = m_numV + num_fixedVertices; //for each fixed vertice, add a new row and col
        MatrixXf m = MatrixXf::Zero(n,n); 
        
		for(int i=0;i<m_numV;i++){ // Eigen::block() not with non constant vars
			for(int j=0;j<m_numV;j++){
				m(i,j) = m_laplaceMatrix(i,j);
			}
		}
        
		// Add 1s in the row and column associated with the fixed point to constain it
        // This will increase L by the size of fixed_verts
        for (int i =0; i< num_fixedVertices; ++i){
			m(m_numV +i, m_fixedVertices[i]) = 1;
			m( m_fixedVertices[i], m_numV +i) = 1;
		}

        m_laplaceMatrix = m;
	}

	void computeDistances(){ //error in here
		m_distances.clear();
        for (int i=0; i< m_numV; i++){
            Vertex v_i = m_vertices[i];
            vector<int> neighbors = getNeighborsOf(i);
            int numNeighbors = neighbors.size();

            vector<Vector4f> distance_i(numNeighbors);

            for (int j=0; j< numNeighbors; j++){
                int nId = neighbors[j];

                Vertex v_j = m_vertices[nId];
                distance_i[j] = Vector4f(v_i.position - v_j.position);
			}
        m_distances.push_back(distance_i);
		}
	}

        
	bool writeMesh(const std::string& filename) {
		// Write off file.
		std::ofstream outFile(filename);
		if (!outFile.is_open()) return false;

		// Write header.
		outFile << "COFF" << std::endl;
		outFile << m_vertices.size() << " " << m_triangles.size() << " 0" << std::endl;

		// Save vertices.
		for (unsigned int i = 0; i < m_vertices.size(); i++) {
			const auto& vertex = m_vertices[i];
			if (vertex.position.allFinite())
				outFile << vertex.position.x() << " " << vertex.position.y() << " " << vertex.position.z() << " "
				<< int(vertex.color.x()) << " " << int(vertex.color.y()) << " " << int(vertex.color.z()) << " " << int(vertex.color.w()) << std::endl;
			else
				outFile << "0.0 0.0 0.0 0 0 0 0" << std::endl;
		}

		// Save faces.
		for (unsigned int i = 0; i < m_triangles.size(); i++) {
			outFile << "3 " << m_triangles[i].idx0 << " " << m_triangles[i].idx1 << " " << m_triangles[i].idx2 << std::endl;
		}

		// Close file.
		outFile.close();

		return true;
	}

	/**
	 * Joins two meshes together by putting them into the common mesh and transforming the vertex positions of
	 * mesh1 with transformation 'pose1to2'. 
	 */
	static SimpleMesh joinMeshes(const SimpleMesh& mesh1, const SimpleMesh& mesh2, Matrix4f pose1to2 = Matrix4f::Identity()) {
		SimpleMesh joinedMesh;
		const auto& vertices1  = mesh1.getVertices();
		const auto& triangles1 = mesh1.getTriangles();
		const auto& vertices2  = mesh2.getVertices();
		const auto& triangles2 = mesh2.getTriangles();

		auto& joinedVertices  = joinedMesh.getVertices();
		auto& joinedTriangles = joinedMesh.getTriangles();

		const unsigned nVertices1 = vertices1.size();
		const unsigned nVertices2 = vertices2.size();
		joinedVertices.reserve(nVertices1 + nVertices2);

		const unsigned nTriangles1 = triangles1.size();
		const unsigned nTriangles2 = triangles2.size();
		joinedTriangles.reserve(nVertices1 + nVertices2);

		// Add all vertices (we need to transform vertices of mesh 1).
		for (int i = 0; i < nVertices1; ++i) {
			const auto& v1 = vertices1[i];
			Vertex v;
			v.position = pose1to2 * v1.position;
			v.color = v1.color;
			joinedVertices.push_back(v);
		}
		for (int i = 0; i < nVertices2; ++i) joinedVertices.push_back(vertices2[i]);

		// Add all faces (the indices of the second mesh need to be added an offset).
		for (int i = 0; i < nTriangles1; ++i) joinedTriangles.push_back(triangles1[i]);
		for (int i = 0; i < nTriangles2; ++i) {
			const auto& t2 = triangles2[i];
			Triangle t{ t2.idx0 + nVertices1, t2.idx1 + nVertices1, t2.idx2 + nVertices1 };
			joinedTriangles.push_back(t);
		}

		return joinedMesh;
	}

	/**
	 * Generates a sphere around the given center point.
	 */
	static SimpleMesh sphere(Vector3f center, float scale = 1.f, Vector4uc color = { 0, 0, 255, 255 }) {
		SimpleMesh mesh;
		Vector4f centerHomogenous = Vector4f{ center.x(), center.y(), center.z(), 1.f };
		
		// These are precomputed values for sphere aproximation.
		const std::vector<double> vertexComponents = { -0.525731, 0, 0.850651 ,0.525731, 0 ,0.850651, -0.525731, 0 ,-0.850651, 0.525731, 0 ,-0.850651, 0, 0.850651, 0.525731, 0, 0.850651, -0.525731, 0, 
			-0.850651, 0.525731, 0, -0.850651, -0.525731, 0.850651, 0.525731, 0, -0.850651, 0.525731, 0, 0.850651, -0.525731, 0, -0.850651, -0.525731, 0 };
		const std::vector<unsigned> faceIndices = { 0, 4, 1, 0, 9, 4, 9, 5, 4, 4, 5, 8, 4, 8, 1, 8, 10, 1, 8, 3, 10, 5, 3, 8, 5, 2, 3, 2, 7, 3, 7, 10,
			3, 7, 6, 10, 7, 11, 6, 11, 0, 6, 0, 1, 6, 6, 1, 10, 9, 0, 11, 9, 11, 2, 9, 2, 5, 7, 2, 11 };

		// Add vertices.
		for (int i = 0; i < 12; ++i) {
			Vertex v;
			v.position = centerHomogenous + scale * Vector4f{ float(vertexComponents[3 * i + 0]), float(vertexComponents[3 * i + 1]), float(vertexComponents[3 * i + 2]), 0.f };
			v.color = color;
			mesh.addVertex(v);
		}

		// Add faces.
		for (int i = 0; i < 20; ++i) {
			mesh.addFace(faceIndices[3 * i + 0], faceIndices[3 * i + 1], faceIndices[3 * i + 2]);
		}

		return mesh;
	}

	/**
	 * Generates a camera object with a given pose.
	 */
	static SimpleMesh camera(const Matrix4f& cameraPose, float scale = 1.f, Vector4uc color = { 255, 0, 0, 255 }) {
		SimpleMesh mesh;
		Matrix4f cameraToWorld = cameraPose.inverse();

		// These are precomputed values for sphere aproximation.
		std::vector<double> vertexComponents = { 25, 25, 0, -50, 50, 100, 49.99986, 49.9922, 99.99993, -24.99998, 25.00426, 0.005185, 
			25.00261, -25.00023, 0.004757, 49.99226, -49.99986, 99.99997, -50, -50, 100, -25.00449, -25.00492, 0.019877 };
		const std::vector<unsigned> faceIndices = { 1, 2, 3, 2, 0, 3, 2, 5, 4, 4, 0, 2, 5, 6, 7, 7, 4, 5, 6, 1, 7, 1, 3, 7, 3, 0, 4, 7, 3, 4, 5, 2, 1, 5, 1, 6 };

		// Add vertices.
		for (int i = 0; i < 8; ++i) {
			Vertex v;
			v.position = cameraToWorld * Vector4f{ scale * float(vertexComponents[3 * i + 0]), scale * float(vertexComponents[3 * i + 1]), scale * float(vertexComponents[3 * i + 2]), 1.f };
			v.color = color;
			mesh.addVertex(v);
		}

		// Add faces.
		for (int i = 0; i < 12; ++i) {
			mesh.addFace(faceIndices[3 * i + 0], faceIndices[3 * i + 1], faceIndices[3 * i + 2]);
		}

		return mesh;
	}

	/**
	 * Generates a cylinder, ranging from point p0 to point p1.
	 */
	static SimpleMesh cylinder(const Vector3f& p0, const Vector3f& p1, float radius, unsigned stacks, unsigned slices, const Vector4uc color = Vector4uc{ 0, 0, 255, 255 }) {
		SimpleMesh mesh;
		auto& vertices = mesh.getVertices();
		auto& triangles = mesh.getTriangles();

		vertices.resize((stacks + 1) * slices);
		triangles.resize(stacks * slices * 2);

		float height = (p1 - p0).norm();

		unsigned vIndex = 0;
		for (unsigned i = 0; i <= stacks; i++)
			for (unsigned i2 = 0; i2 < slices; i2++)
			{
				auto& v = vertices[vIndex++];
				float theta = float(i2) * 2.0f * M_PI / float(slices);
				v.position = Vector4f{ p0.x() + radius * cosf(theta), p0.y() + radius * sinf(theta), p0.z() + height * float(i) / float(stacks), 1.f };
				v.color = color;
			}

		unsigned iIndex = 0;
		for (unsigned i = 0; i < stacks; i++)
			for (unsigned i2 = 0; i2 < slices; i2++) {
				int i2p1 = (i2 + 1) % slices;

				triangles[iIndex].idx0 = (i + 1) * slices + i2;
				triangles[iIndex].idx1 = i * slices + i2;
				triangles[iIndex].idx2 = i * slices + i2p1;

				triangles[iIndex + 1].idx0 = (i + 1) * slices + i2;
				triangles[iIndex + 1].idx1 = i * slices + i2p1;
				triangles[iIndex + 1].idx2 = (i + 1) * slices + i2p1;

				iIndex += 2;
			}

		Matrix4f transformation = Matrix4f::Identity();
		transformation.block(0, 0, 3, 3) = face(Vector3f{ 0, 0, 1 }, p1 - p0);
		transformation.block(0, 3, 3, 1) = p0;
		mesh.transform(transformation);

		return mesh;
	}

private:
	vector<Vertex> m_vertices;
	vector<int> m_fixedVertices;
	vector<Vertex> m_verticesPrime;
	vector<Triangle> m_triangles;
	MatrixXf m_neighborMatrix;
	vector<MatrixXf> m_cellRotations;
	MatrixXf m_edgeMatrix;
	MatrixXf m_laplaceMatrix;
	vector<std::vector<unsigned int>> m_verticesToFaces;
	MatrixXf m_weightMatrix, m_weightSum;
	vector<vector<Vector4f>> m_distances;
	int m_numV;

	/**
	 * Returns a rotation that transforms vector vA into vector vB.
	 */
	static Matrix3f face(const Vector3f& vA, const Vector3f& vB) {
		auto a = vA.normalized();
		auto b = vB.normalized();
		auto axis = b.cross(a);
		float angle = acosf(a.dot(b));
		
		if (angle == 0.0f) {  // No rotation
			return Matrix3f::Identity();
		}

		// Convert the rotation from SO3 to matrix notation.
		// First we create a skew symetric matrix from the axis vector.
		Matrix3f skewSymetricMatrix;
		skewSymetricMatrix.setIdentity();
		skewSymetricMatrix(0, 0) = 0;			skewSymetricMatrix(0, 1) = -axis.z();	skewSymetricMatrix(0, 2) = axis.y();
		skewSymetricMatrix(1, 0) = axis.z();	skewSymetricMatrix(1, 1) = 0;			skewSymetricMatrix(1, 2) = -axis.x();
		skewSymetricMatrix(2, 0) = -axis.y();	skewSymetricMatrix(2, 1) = axis.x();	skewSymetricMatrix(2, 2) = 0;

		// We compute a rotation matrix using Rodrigues formula.
		Matrix3f rotation = Matrix3f::Identity() + sinf(angle) * skewSymetricMatrix + (1 - cos(angle)) * skewSymetricMatrix * skewSymetricMatrix;

		return rotation;
	}
};
