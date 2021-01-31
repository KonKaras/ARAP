#pragma once

#include <iostream>
#include <fstream>

#include "Eigen.h"
using namespace std;

struct Vertex {
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


class SimpleMesh {
public:
	SimpleMesh() {}

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

	bool loadMesh(const std::string& filename) {
		// Read off file (Important: Only .off files are supported).
		m_vertices.clear();
		m_vertices_prime.clear();
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
		unsigned int num_v = 0; //vertices
		unsigned int num_p = 0; //faces
		unsigned int num_e = 0; //edges
		file >> num_v >> num_p >> num_e;

		m_num_v = num_v;
		m_num_p = num_p;

		// Read vertices.
		if (std::string(string1).compare("COFF") == 0) {
			// We have color information.
			for (unsigned int i = 0; i < num_v; i++) {
				Vertex v;
				file >> v.position.x() >> v.position.y() >> v.position.z();

				// Colors are stored as integers. We need to convert them.
				Vector4i colorInt;
				file >> colorInt.x() >> colorInt.y() >> colorInt.z() >> colorInt.w();
				v.color = Vector4uc((unsigned char)colorInt.x(), (unsigned char)colorInt.y(), (unsigned char)colorInt.z(), (unsigned char)colorInt.w());
				m_vertices.push_back(v);
				m_vertices_prime.push_back(v);
			}
		}
		else if (std::string(string1).compare("OFF") == 0) {
			// We only have vertex information.
			for (unsigned int i = 0; i < num_v; i++) {
				Vertex v;
				file >> v.position.x() >> v.position.y() >> v.position.z();
				v.color.x() = 0;
				v.color.y() = 0;
				v.color.z() = 0;
				v.color.w() = 255;
				m_vertices.push_back(v);
				m_vertices_prime.push_back(v);
			}
		}
		else {
			std::cout << "Incorrect mesh file type." << std::endl;
			return false;
		}

		m_neighbor_matrix = MatrixXf::Zero(m_num_v, m_num_v); 
		m_vertices_to_faces = std::vector<std::vector<unsigned int>>(m_num_v); //vtf[i] contains the face IDs that contain vertex ID i
		m_triangles = vector<Triangle>(num_p);

		// Read faces (i.e. triangles).
		for (unsigned int i = 0; i < num_p; i++) {
			unsigned int num_vs;
			file >> num_vs;
			ASSERT(num_vs == 3 && "We can only read triangular mesh.");

			Triangle t;
			file >> t.idx0 >> t.idx1 >> t.idx2;
			m_triangles[i] = t;

			addFaceToAdjacencyMatrix(t.idx0, t.idx1, t.idx2);  // 1 for neighbors, 0 else
			(m_vertices_to_faces[t.idx0]).push_back(i); // list of facenumbers for each vertex
			(m_vertices_to_faces[t.idx1]).push_back(i);
			(m_vertices_to_faces[t.idx2]).push_back(i);

		}

		cout << "Mesh has "<<num_v<<" vertices, "<<num_p<<" faces and "<<num_e<<" edges"<<endl;
		return true;
	}

	void printPs() {
		cout << "m_vertices:" << endl;
		for (Vertex v : m_vertices) {
			cout << v.position.x() << " " << v.position.y() << " " << v.position.z() << endl;
		}
	}

	void printPPrimes() {
		cout << "m_vertices_prime:" << endl;
		for (Vertex v : m_vertices_prime) {
			cout << v.position.x() << " " << v.position.y() << " " << v.position.z() << endl;
		}
	}

	void addFaceToAdjacencyMatrix(int v1, int v2, int v3) {
		m_neighbor_matrix(v1, v2) = 1;
		m_neighbor_matrix(v2, v1) = 1;
		m_neighbor_matrix(v1, v3) = 1;
		m_neighbor_matrix(v3, v1) = 1;
		m_neighbor_matrix(v2, v3) = 1;
		m_neighbor_matrix(v3, v2) = 1;
		m_nonzero_neighbor_matrix += 6;
	}

	int getNonZeroNeighborMatrixEntries(){
		return m_nonzero_neighbor_matrix;
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
		for (int i = 0; i < m_num_v; i++) {
			if (m_neighbor_matrix(vertId, i) == 1)
				neighbors.push_back(i);
		}
		return neighbors;
	}

	Vector3f getVertex(int i) {
		return m_vertices[i].position;
	}

	Vector3f getDeformedVertex(int i) {
		return m_vertices_prime[i].position;
	}

	int getNumberOfVertices() {
		return m_num_v;
	}

	int getNumberOfFaces() {
		return m_num_p;
	}

	void setPPrime(MatrixXf pprime) {
		for (int i = 0; i < m_num_v; i++) {
			Vertex v;
			v.position.x() = pprime(i, 0);
			v.position.y() = pprime(i, 1);
			v.position.z() = pprime(i, 2);
			m_vertices_prime[i] = v;
		}
	}

	void setPPrime(int index, Vector3f position) {
		m_vertices_prime[index].position = position;
	}


	void copyPPrime() {
		for (int i = 0; i < m_num_v; i++) {
			m_vertices[i] = m_vertices_prime[i];
		}
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

	float computeUniformWeightForVertex(int i){
		vector<int> neighbors_i = getNeighborsOf(i);
		float wij = 1 / (float) neighbors_i.size();
		cout<<"Uniform weight for vertex "<<i<<" with "<<neighbors_i.size() <<" is : "<<wij<<endl;
		return wij;
	}

	float computeCotangentWeightForPair(int i, int j) {
		vector<Triangle> local_faces;
		for (unsigned int id : m_vertices_to_faces[i]) {
			Triangle f = m_triangles[id];
			if (partOfFace(i, j, f)) 
				local_faces.push_back(f);
		}
		assert(local_faces.size() <= 2);

		Vertex vertex_i = m_vertices[i];
		Vertex vertex_j = m_vertices[j];

		// Using the following weights: 0.5 * (cot(alpha) + cot(beta))
		float cot_theta_sum = 0.0f;

		for (Triangle f : local_faces) {
			int other_vertex_id = getThirdFacePoint(i, j, f);
			Vertex vertex_o = m_vertices[other_vertex_id];
			float theta = cotan(vertex_i.position - vertex_o.position, vertex_j.position - vertex_o.position); 
			// float theta = angleBetweenVectors(vertex_i.position - vertex_o.position, vertex_j.position - vertex_o.position); 
			cout<<"Theta for i "<<i<<" and j "<<j<<" and other v "<<other_vertex_id<<" = "<<theta<<endl;
			// cot_theta_sum += (1 / tan(theta));
			cot_theta_sum += theta;
		}

		return cot_theta_sum * 0.5f;
	}

	float cotan(Vector3f a, Vector3f b){
    	return (a.dot(b)) / (a.cross(b)).norm();
	}

	double angleBetweenVectors(Vector3f a, Vector3f b) {
		return acos((a).dot(b) / (a.norm() * b.norm())) * 180 / M_PI;
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
	vector<Vertex> m_vertices_prime;
	vector<Triangle> m_triangles;
	MatrixXf m_neighbor_matrix;
	vector<std::vector<unsigned int>> m_vertices_to_faces;
	int m_num_v;
	int m_num_p;
	int m_nonzero_neighbor_matrix = 0;


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
