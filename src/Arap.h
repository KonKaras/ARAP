#include <iostream>
#include <fstream>

#include "Eigen.h"
#include "SimpleMesh.h"
using namespace std;

void estimateRotation(SimpleMesh mesh, int vertexID) {
		// Assume vertices are fix, solve for rotations with Procrustes

        // Sorkine paper: PDP' for vertex i with neighbors j (P contains edges of fan as cols, P' edges of deformed fan, D diagonal with weights w_ij)
		Matrix3f rotation = Matrix3f::Identity();
        vector<int> neighbors = mesh.getNeighborsOf(vertexID);
        int numNeighbors = neighbors.size();

		MatrixXf P (3, numNeighbors);
		MatrixXf PPrime (3, numNeighbors);
        MatrixXf D (numNeighbors, numNeighbors);

        cout<< "In estimate rots()"<< endl;

		for ( int j = 0; j< numNeighbors; ++j)
		{   
            int neighborVertex = neighbors[j];
			P.col(j) = mesh.getVertex(vertexID) - mesh.getVertex(neighborVertex);
			PPrime.col(j) = mesh.getDeformedVertex(vertexID) - mesh.getDeformedVertex(neighborVertex);
            D(j,j) = mesh.getWeight(vertexID, neighborVertex);
		} 

        cout<< "ESTIMATEROTS "<<vertexID<<endl;
        cout<< P << endl;
        cout<< " "<<endl;
        cout<< PPrime << endl;
        cout<< " "<<endl;
        cout<< D << endl;
        cout<< " "<<endl;

		JacobiSVD<MatrixXf> svd(P * D * PPrime.transpose(), ComputeFullU | ComputeFullV);

		rotation = svd.matrixV() * svd.matrixU().transpose();

		if(rotation.determinant() == -1)
		{
            MatrixXf svd_u = svd.matrixU();
            svd_u.rightCols(1) = svd_u.rightCols(1) * -1; 
			//rotation = svd.matrixV() * Matrix3f::Identity() * Vector3f(1,1,-1) * svd.matrixU().transpose();
			rotation = svd.matrixV() * svd_u.transpose();
		}
        cout<< rotation << endl;

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
    vector<int> fixedVertices = mesh.getFixedVertices();
    for(int i=0; i< mesh.getNumberOfFixedVertices(); ++i){
        int fixedVertex = fixedVertices[i];
        b.row(mesh.getNumberOfVertices()+i) = mesh.getVertexForFillingB(fixedVertex);
    }

    cout<<"b: \n"<<b<<endl;
    //Solve LES with Cholesky, L positive definite // TODO test sparse cholesky on sparse eigen matrices
    cout<<"Solving LES ..." <<endl;
    MatrixXf PPrime = mesh.getLaplaceMatrix().llt().solve(b);
    cout<<"Done!"<<endl;
    cout<<"pprime: \n"<<PPrime<<endl;
    mesh.setPPrime(PPrime);
}

float calculateEnergy(SimpleMesh mesh){ // TODO: not sure if implemented energy function correctly
    float energy= 0.0f;
    for(int i=0; i<mesh.getNumberOfVertices(); ++i){
        vector<int> neighbors = mesh.getNeighborsOf(i);
        // int numNeighbors = neighbors.size();
        float cell_energy=0;

        for ( int j : neighbors)
        {
            Vector3f v= ((mesh.getDeformedVertex(i) - mesh.getDeformedVertex(j)) - mesh.getRotation(i) * (mesh.getVertex(i) - mesh.getVertex(j))); 
            cell_energy += mesh.getWeight(i,j) * pow(v[0]*v[0] + v[1]*v[1]+v[2]*v[2], 2);
        } 
        energy+=cell_energy;
    }
    return energy;
}

void applyDeformation(SimpleMesh mesh, int handleID, Vector4f handleNewPosition, int iterations){
    mesh.setNewHandlePosition(handleNewPosition);
    float energy=0.0f;
    int iter=0;
    cout<<"Applying deformation for handle with ID " << handleID << " to new position " << handleNewPosition.x() <<","<< handleNewPosition.y()<< ","<< handleNewPosition.z()<<endl;
    while(iter<iterations){
        cout<<"[Iteration "<<iter<<"]"<<endl;

        //TODO Initial guess for pprime ?

        for(int i=0; i< mesh.getNumberOfVertices(); ++i){
            estimateRotation(mesh, i);
        }
        estimateVertices(mesh);
        float energy_i = calculateEnergy(mesh);        
        cout<< "Iteration: "<< iter<< "  Local error: "<< energy_i << endl;

        mesh.copyPPrime();

        energy = energy_i;
        iter++;
    }
    cout << "Resulting energy: "<< energy<< endl;
    cout << "PPrime[handleID] is "<< mesh.getDeformedVertex(handleID).x() <<","<< mesh.getDeformedVertex(handleID).y()<< ","<< mesh.getDeformedVertex(handleID).z()<<endl;
    // assert(mesh.getDeformedVertex(handleID).x() == handleNewPosition.x());
    // assert(mesh.getDeformedVertex(handleID).y() == handleNewPosition.y());
    // assert(mesh.getDeformedVertex(handleID).z() == handleNewPosition.z());
}