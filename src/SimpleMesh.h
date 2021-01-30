#pragma once

#include <iostream>
#include <fstream>

#include "Eigen.h"
using namespace std;

struct Vertex {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		// Position stored as 3 floats
		Vector3f position;
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


struct Constraint {
	int vertexID;
	Vector3f position;
};


class SimpleMesh {
public:
	SimpleMesh(){
	}
	MatrixXf m_b;

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

	bool loadMeshFromGUI(){
		return true;
	}

	bool loadMesh(const std::string& filename, vector<int> fixedPoints) {
		// Read off file (Important: Only .off files are supported).
		m_vertices.clear();
		m_verticesPrime.clear();
		m_triangles.clear();
		m_constraints.clear();
		m_verticesToFaces.clear();
		m_cellRotations.clear();
		m_precomputedPMatrices.clear();
		
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

		m_numV = numV;

		// // Fixed vertices
		//  for(int f : fixedPoints){
		//  	m_fixedVertices.push_back(f);
		//  }

		//Handle is last fixed Vertex. Handle darf sich auch nicht bewegen, da er ja eine fixe zielposition zugewiesen bekommen hat
		//m_handleID = handle;
		
		// Read vertices.
		if (std::string(string1).compare("COFF") == 0) {
			// We have color information.
			for (unsigned int i = 0; i < numV; i++) {
				Vertex v;
				file >> v.position.x() >> v.position.y() >> v.position.z();

				// Colors are stored as integers. We need to convert them.
				Vector4i colorInt;
				file >> colorInt.x() >> colorInt.y() >> colorInt.z() >> colorInt.w();
				v.color = Vector4uc((unsigned char)colorInt.x(), (unsigned char)colorInt.y(), (unsigned char)colorInt.z(), (unsigned char)colorInt.w());
				m_vertices.push_back(v);
				m_verticesPrime.push_back(v);
			}

			for (int i : fixedPoints) {
				Constraint c;
				c.vertexID = i;
				c.position = m_vertices[i].position;
				m_constraints.push_back(c);
				// cout<<"Added constraint with ID "<<i<<" and position "<<v.position<<" Constraints length = "<<m_constraints.size()<<endl;
			}
		}
		else if (std::string(string1).compare("OFF") == 0) {
			// We only have vertex information.
			for (unsigned int i = 0; i < numV; i++) {
				Vertex v;
				file >> v.position.x() >> v.position.y() >> v.position.z();

				v.color.x() = 0;
				v.color.y() = 0;
				v.color.z() = 0;
				v.color.w() = 255;

				m_vertices.push_back(v);
				m_verticesPrime.push_back(v);

			}

			for (int i : fixedPoints) {
				//if (i >= 0) {
					Constraint c;
					c.vertexID = i;
					c.position = m_vertices[i].position;
					m_constraints.push_back(c);
				//}
				//cout << "Added constraint with ID " << i << " and position " << m_vertices[i].position << " Constraints length = " << m_constraints.size() << endl;
			}

			for (Constraint i : m_constraints) {
				cout << i.vertexID << endl;
			}
		}
		else {
			std::cout << "Incorrect mesh file type." << std::endl;
			return false;
		}
		

		//Speicherallokation der wichtigen Matrizen und Listen
		m_neighborMatrix = MatrixXf::Zero(m_numV, m_numV); // Element (i,j)=1, wenn i und j nachbarb, sonst 0
		m_verticesToFaces = std::vector<std::vector<unsigned int>>(m_numV); //vtf[i] contains the face IDs that contain vertex ID i
		m_cellRotations = std::vector<MatrixXf>(m_numV); // eine rotationsmatrix pro vertex
		m_precomputedPMatrices = vector<MatrixXf>(m_numV); // Vorberechnung der P matrix (paper seite 4 anfang)
		m_triangles = vector<Triangle>(numP);
		m_b = MatrixXf::Zero(m_numV, 3);

		// Read faces (i.e. triangles).
		for (unsigned int i = 0; i < numP; i++) {
			unsigned int num_vs;
			file >> num_vs;
			ASSERT(num_vs == 3 && "We can only read triangular mesh.");

			Triangle t;
			file >> t.idx0 >> t.idx1 >> t.idx2;
			m_triangles[i] = t;

			//fill adjacency matrix
			addFaceToAdjacencyMatrix(t.idx0, t.idx1, t.idx2);  // 1 for neighbors, 0 else
			(m_verticesToFaces[t.idx0]).push_back(i); // list of facenumbers for each vertex
			(m_verticesToFaces[t.idx1]).push_back(i);
			(m_verticesToFaces[t.idx2]).push_back(i);

		}

		for (int i = 0; i < numV; i++) {
			m_cellRotations[i] = MatrixXf::Zero(3, 3);
		}

		cout << "numVertices: " << numV << endl;
		cout << "numFaces: " << numP << endl;
		cout << "numEdges: " << numE << endl;

		// cout << "Adjacencymatrix: " << m_neighborMatrix <<endl;
		// cout << m_neighborMatrix.rows() << " - " << m_neighborMatrix.cols()<< endl;

		buildWeightMatrix(); // Die weights werden hier vorberechnet

		// cout << "m_weightmatrix: " << m_weightMatrix <<endl;
		// cout << "m_weightSum: " << m_weightSum <<endl;
		// computeDistances();
		calculateSystemMatrix(); // Die Matrix L (paper seite 5 anfang) wird hier berechnet, also die linke seite des LGS. Siehe auch gegebenen ARAP code.
		precomputePMatrix(); // Die P matrix (seite 4 anfang) wird berechnet

		return true;
	}

	// void applyConstrainedPoints(int handleID, vector<int> fixedPoints){
	// 	m_handleID = handleID;
	// 	cout << "fixedpoints given to apply " << fixedPoints.size()<<endl;
	// 	// Fixed vertices
	// 	for(int f : fixedPoints){
	// 		m_fixedVertices.push_back(f);
	// 	}
	// }

	void setHandleConstraint(int handleID, Vector3f newHandlePosition) {
		for (int i = 0; i < m_constraints.size(); ++i) {
			if (m_constraints[i].vertexID == handleID) {
				m_constraints[i].position = newHandlePosition;
				//cout << "Set fixed vertex with ID " << m_constraints[i].vertexID << " to " << m_constraints[i].position << endl;
			}
		}
	}


	Vector3f getConstraintI(int id) {
		for (Constraint c : m_constraints) {
			if (c.vertexID == id) {
				return c.position;
			}
		}
		throw std::invalid_argument("received id which is not in fixed vertices");
	}

	void printPs() {
		cout << "m_vertices:" << endl;
		for (Vertex v : m_vertices) {
			cout << v.position.x() << " " << v.position.y() << " " << v.position.z() << endl;
		}
	}

	void printNewHandlePosition() {
		cout << "New handle position: " << m_newHandlePosition.x() << " " << m_newHandlePosition.y() << " " << m_newHandlePosition.z() << endl;
	}

	void printPPrimes() {
		cout << "m_verticesPrime:" << endl;
		for (Vertex v : m_verticesPrime) {
			cout << v.position.x() << " " << v.position.y() << " " << v.position.z() << endl;
		}
	}

	//Precompute P_i s for the arap estimateRotations() function (page 4 beginning)
	void precomputePMatrix() {

		for (int i = 0; i < m_numV; i++) {
			vector<int> neighbors = getNeighborsOf(i);
			int numNeighbors = neighbors.size();
			MatrixXf P = MatrixXf::Zero(3, numNeighbors);
			for (int j = 0; j < numNeighbors; j++)
			{
				int neighborVertex = neighbors[j];
				P.col(j) = getVertex(i) - getVertex(neighborVertex);
			}
			m_precomputedPMatrices[i] = P;
		}
	}

	MatrixXf getPrecomputedP(int i) {
		return m_precomputedPMatrices[i];
	}

	void addFaceToAdjacencyMatrix(int v1, int v2, int v3) {
		m_neighborMatrix(v1, v2) = 1;
		m_neighborMatrix(v2, v1) = 1;
		m_neighborMatrix(v1, v3) = 1;
		m_neighborMatrix(v3, v1) = 1;
		m_neighborMatrix(v2, v3) = 1;
		m_neighborMatrix(v3, v2) = 1;
	}

	//retrun true if both i and j are part of face f
	bool partOfFace(int i, int j, Triangle f) {
		return (f.idx0 == i && (f.idx1 == j || f.idx2 == j)) ||
			(f.idx1 == i && (f.idx0 == j || f.idx2 == j)) ||
			(f.idx2 == i && (f.idx1 == j || f.idx0 == j));
	}

	//return list of all neighbors of vertex i
	vector<int> getNeighborsOf(int vertId) {
		vector<int> neighbors;
		for (int i = 0; i < m_numV; i++) {
			if (m_neighborMatrix(vertId, i) == 1)
				neighbors.push_back(i);
		}
		return neighbors;
	}

	Vector3f getVertex(int i) {
		return m_vertices[i].position;
	}

	// When filling the B vector the new handle position is needed, however for other vertices their original position is needed which is 
	// Why we have the two functions getVertex() and getVertexForFillingB(), which differs between handle and not handle
	//TODO might not be correct to fill row of handle in b vector this way...
	Vector3f getVertexForFillingB(int i) {
		if (i == m_handleID)
			return m_newHandlePosition;
		else
			return m_vertices[i].position;
	}

	Vector3f getDeformedVertex(int i) {
		return m_verticesPrime[i].position;
	}

	float getWeight(int i, int j) {
		return m_weightMatrix(i, j);
	}

	int getNumberOfVertices() {
		return m_numV;
	}

	int getNumberOfConstraints() {
		return m_constraints.size();
	}

	MatrixXf getRotation(int i) {
		return m_cellRotations[i];
	}

	void setRotation(int i, MatrixXf r) {
		m_cellRotations[i] = r;
	}

	void setPPrime(MatrixXf pprime) {
		for (int i = 0; i < m_numV; i++) {
			Vertex v;
			v.position.x() = pprime(i, 0);
			v.position.y() = pprime(i, 1);
			v.position.z() = pprime(i, 2);
			m_verticesPrime[i] = v;
		}
	}

	void setPPrime(int index, Vector3f position) {
		m_verticesPrime[index].position = position;
	}


	void copyPPrime() {
		for (int i = 0; i < m_numV; i++) {
			m_vertices[i] = m_verticesPrime[i];
		}
	}

	MatrixXf getSystemMatrix() {
		return m_systemMatrix;
	}

	int getThirdFacePoint(int i, int j, Triangle f) {
		vector<int> points;
		points.push_back(f.idx0);
		points.push_back(f.idx1);
		points.push_back(f.idx2);
		for (int idx : points) {
			if (idx != i && idx != j) {
				return idx;
			}
		}
		return -1;
	}

	void setNewHandlePosition(Vector3f newPosition) {
		m_newHandlePosition = newPosition;
	}

	void buildWeightMatrix() {
		cout << "Generating Weight Matrix" << endl;
		//TODO compute weights
		m_weightMatrix = MatrixXf::Ones(m_numV, m_numV);
		//m_weightMatrix = MatrixXf::Zeros(m_numV, m_numV);
		//m_weightSum = MatrixXf::Zero(m_numV, m_numV);
		/*
		for (int vertex_id = 0; vertex_id < m_numV; vertex_id++) {
			vector<int> neighbors = getNeighborsOf(vertex_id);
			for (int neighbor_id : neighbors)
				assignWeightForPair(vertex_id, neighbor_id);
		}
		*/
		
	}

	void assignWeightForPair(int i, int j) {
		float weightIJ;
		if (m_weightMatrix(j, i) == 0) //If the opposite weight has not been computed, do so
			weightIJ = computeWeightForPair(i, j);
		else
			weightIJ = m_weightMatrix(j, i);

		//m_weightSum(i, i) += (weightIJ * 0.5);
		//m_weightSum(j, j) += (weightIJ * 0.5);
		m_weightMatrix(i, j) = weightIJ;
	}

	double angleBetweenVectors(Vector3f a, Vector3f b) {
		return acos((a).dot(b) / (a.norm() * b.norm())) * 180 / M_PI;
	}

	vector<Constraint> getConstraints() {
		return m_constraints;
	}

	bool isInConstraints(int i) {
		for (Constraint c : m_constraints) {
			if (c.vertexID == i) return true;
		}
		return false;
	}

	bool isHandle(int i) {
		return i == m_handleID;
	}

	Vector3f getHandleNewPosition() {
		return m_newHandlePosition;
	}

	float computeWeightForPair(int i, int j) {
		vector<Triangle> localFaces;
		for (unsigned int id : m_verticesToFaces[i]) {
			Triangle f = m_triangles[id];
			if (partOfFace(i, j, f)) //If the face contains both I and J, add it
				localFaces.push_back(f);
		}
		// Either a normal face or a boundry edge, otherwise bad mesh
		assert(localFaces.size() <= 2);

		Vertex vertex_i = m_vertices[i];
		Vertex vertex_j = m_vertices[j];

		// Using the following weights: 0.5 * (cot(alpha) + cot(beta))

		float cot_theta_sum = 0.0f;

		for (Triangle f : localFaces) {
			int other_vertex_id = getThirdFacePoint(i, j, f);
			Vertex vertex_o = m_vertices[other_vertex_id];

			float theta = angleBetweenVectors(vertex_i.position - vertex_o.position, vertex_j.position - vertex_o.position); //TODO is this correct???
			// cout<<"Theta for i "<<i<<" and j "<<j<<" and other v "<<other_vertex_id<<" = "<<theta<<endl;
			cot_theta_sum += (1 / tan(theta));
		}

		return cot_theta_sum * 0.5f;
	}


	void calculateSystemMatrix() {

		m_systemMatrix = MatrixXf::Zero(m_numV, m_numV);
		for (int i = 0; i < m_numV; i++) {
			vector<int> neighbors = getNeighborsOf(i);
			int numNeighbors = neighbors.size();
			for (int j = 0; j < numNeighbors; ++j)
			{
				int neighborVertex = neighbors[j];
				m_systemMatrix(i, i) += m_weightMatrix(i, j);
				m_systemMatrix(i, neighborVertex) = -m_weightMatrix(i, neighborVertex);
			}
		}

		for (Constraint c : m_constraints) {
			int i = c.vertexID;
			m_systemMatrix.row(i).setZero();
			m_systemMatrix(i, i) = 1;
		}
	}

	// Writes mesh to file
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

private:
	vector<Vertex> m_vertices;
	// vector<int> m_fixedVertices;
	// vector<Vertex> m_fixedVerticesPositions;
	vector<Vertex> m_verticesPrime;
	vector<Triangle> m_triangles;
	MatrixXf m_neighborMatrix;
	vector<MatrixXf> m_cellRotations;
	MatrixXf m_systemMatrix;
	vector<std::vector<unsigned int>> m_verticesToFaces;
	MatrixXf m_weightMatrix;//, m_weightSum;
	// vector<vector<Vector4f>> m_distances;
	int m_numV;
	int m_handleID;
	Vector3f m_newHandlePosition;
	vector<MatrixXf> m_precomputedPMatrices;
	vector<Constraint> m_constraints;


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
