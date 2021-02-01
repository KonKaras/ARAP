#pragma once

#include <iostream>
#include <fstream>
#include <Eigen/Sparse>

#include "Eigen.h"
#include "SimpleMesh.h"
using namespace std;

#define THRESHOLD 1.0e-3f

struct Constraint {
	int vertexID;
	Vector3d position;
};

class ArapDeformer {

public:
    ArapDeformer();
	ArapDeformer(SimpleMesh *mesh);
    void initDeformation(vector<int> fixed_points);
    void setHandleConstraint(int handleID, Vector3d newHandlePosition);
    void applyDeformation(vector<int> fixed_points, int handleID, Vector3d handleNewPosition, int iterations);
    SimpleMesh m_mesh;

private:

    void estimateRotation();
    void updateB();
    void estimateVertices();
    double calculateEnergy();
    void buildWeightMatrix(); 
    void calculateSystemMatrix();
    Vector3d getConstraintI(int id);
    bool isInConstraints(int i);
    void updateSystemMatrix();

    //SimpleMesh m_mesh;
    vector<MatrixXd> m_cell_rotations;
    MatrixXd m_system_matrix;
    MatrixXd m_system_matrix_original;
    MatrixXd m_weight_matrix;
    SparseMatrix<double> m_system_matrix_sparse;
    int m_num_v;
    int m_num_p;
    MatrixXd m_b;
    int m_handle_id;
    Vector3d m_new_handle_position;
    vector<Constraint> m_constraints;

};

