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
	Vector3f position;
};

class ArapDeformer {

public:
    ArapDeformer();
	ArapDeformer(SimpleMesh *mesh);
    void initDeformation(vector<int> fixed_points);
    void setHandleConstraint(int handleID, Vector3f newHandlePosition);
    void applyDeformation(int handleID, Vector3f handleNewPosition, int iterations);

private:

    void estimateRotation();
    void updateB();
    void estimateVertices();
    float calculateEnergy();
    void buildWeightMatrix(); 
    void calculateSystemMatrix();
    Vector3f getConstraintI(int id);
    bool isInConstraints(int i);

    SimpleMesh m_mesh;
    vector<MatrixXf> m_cell_rotations;
    MatrixXf m_system_matrix;
    MatrixXf m_weight_matrix;
    int m_num_v;
    int m_num_p;
    MatrixXf m_b;
    int m_handle_id;
    Vector3f m_new_handle_position;
    vector<Constraint> m_constraints;

};

