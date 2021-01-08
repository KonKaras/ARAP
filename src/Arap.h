#include <iostream>
#include <fstream>

#include "Eigen.h"
#include "SimpleMesh.h"
using namespace std;

void estimateRotation(SimpleMesh mesh, int vertexID) {
		// Assume vertices are fix, solve for rotations

        // Sorkine paper: PDP' for vertex i with neighbors j (P contains edges of fan as cols, P' edges of deformed fan, D diagonal with weights w_ij)
		Matrix3f rotation = Matrix3f::Identity();
        vector<int> neighbors = mesh.getNeighborsOf(vertexID);
        int numNeighbors = neighbors.size();

		MatrixXf P (3, numNeighbors);
		MatrixXf PPrime (3, numNeighbors);
        MatrixXf D (numNeighbors, numNeighbors);

		for ( int j = 0; j< numNeighbors; ++j)
		{
			P.col(j) = mesh.getVertex(vertexID) - mesh.getVertex(j);
			PPrime.col(j) = mesh.getDeformedVertex(vertexID) - mesh.getDeformedVertex(j);
            D(j,j) = mesh.getWeight(vertexID, j);
		} 
		JacobiSVD<MatrixXf> svd(P * D * PPrime.transpose(), ComputeFullU | ComputeFullV);
		rotation = svd.matrixV() * svd.matrixU().transpose();

		if(rotation.determinant() == -1)
		{
            MatrixXf svd_u = svd.matrixU();
            svd_u.rightCols(1) = svd_u.rightCols(1) * -1;
			//rotation = svd.matrixV() * Matrix3f::Identity() * Vector3f(1,1,-1) * svd.matrixU().transpose();
			rotation = svd.matrixV() * svd_u.transpose();
		}

		mesh.setRotation(vertexID, rotation);
}

Vector3f calculateB(SimpleMesh mesh, int i){
    vector<int> neighbors = mesh.getNeighborsOf(i);
    int numNeighbors = neighbors.size();
    Vector3f b = Vector3f::Zero();
    for ( int j = 0; j< numNeighbors; ++j)
    {
        b += mesh.getWeight(i, j) * 0.5 * (mesh.getRotation(i)+ mesh.getRotation(j))*(mesh.getVertex(i) - mesh.getVertex(j));
    }
    return b;
}

void estimateVertices(SimpleMesh mesh){
    MatrixXf b(mesh.getNumberOfVertices()+ mesh.getNumberOfFixedVertices(), 3);
    for ( int i = 0; i< mesh.getNumberOfVertices(); ++i)
    {
        b.row(i) = calculateB(mesh, i);
    }
    for(int i=0; i< mesh.getNumberOfFixedVertices(); ++i){
        b.row(mesh.getNumberOfVertices()+i) = mesh.getVertex(i);
    }

    //Solve LES with Cholesky, L positive definite // TODO test sparse cholesky on sparse eigen matrices
    MatrixXf PPrime = mesh.getLaplaceMatrix().llt().solve(b);
    mesh.setPPrime(PPrime);
}

float calculateEnergy(SimpleMesh mesh){ // TODO: not sure if implemented energy function correctly
    cout<<"calculateenergy"<<endl;
    float energy= 0.0f;
    for(int i=0; i<mesh.getNumberOfVertices(); ++i){
        vector<int> neighbors = mesh.getNeighborsOf(i);
        int numNeighbors = neighbors.size();
        float energy=0;

        for ( int j = 0; j< numNeighbors; ++j)
        {
            Vector3f v= ((mesh.getVertex(i) - mesh.getVertex(j)) - mesh.getRotation(i) * (mesh.getDeformedVertex(i) - mesh.getDeformedVertex(j))); 
            energy += pow(v[0]*v[0] + v[1]*v[1]+v[2]*v[2], 2);
        } 
    }
    return energy;
}

void applyDeformation(SimpleMesh mesh, int iterations){
    float energy=0.0f;
    cout<<"Applying deformation"<<endl;
    while(iterations>0){

        //TODO Initial guess for pprime needed

        for(int i=0; i< mesh.getNumberOfVertices(); ++i){
            estimateRotation(mesh, i);
        }
        estimateVertices(mesh);
        float energy_i = calculateEnergy(mesh);        
        cout<< "Iterations left: "<< iterations<< "  Local error: "<< energy_i << endl;

        energy = energy_i;
        iterations--;
    }
    cout<< "Resulting energy: "<< energy<< endl;
}