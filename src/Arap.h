#include <iostream>
#include <fstream>

#include "Eigen.h"
#include "SimpleMesh.h"
using namespace std;

void estimateRotation(SimpleMesh *mesh, int vertexID) {
		// Assume vertices are fix, solve for rotations with Procrustes

        // Sorkine paper: PDP' for vertex i with neighbors j (P contains edges of fan as cols, P' edges of deformed fan, D diagonal with weights w_ij)
		Matrix3f rotation = Matrix3f::Identity();
        vector<int> neighbors = mesh->getNeighborsOf(vertexID);
        int numNeighbors = neighbors.size();

		// P precomputed in SimpleMesh.loadMesh()
		MatrixXf PPrime =MatrixXf::Zero(3, numNeighbors); // TODO did I miss any initialization of the deformed vertices?
        MatrixXf D = MatrixXf::Zero(numNeighbors, numNeighbors);

		for ( int j = 0; j< numNeighbors; ++j)
		{   
            int neighborVertex = neighbors[j];
			PPrime.col(j) = mesh->getDeformedVertex(vertexID) - mesh->getDeformedVertex(neighborVertex);
            D(j,j) = mesh->getWeight(vertexID, neighborVertex);
		} 

        MatrixXf P = mesh->getPrecomputedP(vertexID);

        // cout<< "ESTIMATEROTS "<<vertexID<<endl;
        // cout<< P << endl;
        // cout<< " "<<endl;
        // cout<< PPrime << endl;
        // cout<< " "<<endl;
        // cout<< D << endl;
        // cout<< " "<<endl;

		JacobiSVD<MatrixXf> svd(P * D * PPrime.transpose(), ComputeThinU | ComputeThinV);

		rotation = svd.matrixV().transpose() * svd.matrixU(); // TODO not sure which, U or V or both, need to be transposed

		if(rotation.determinant() < 0)
		{
            MatrixXf svd_u = svd.matrixU();
            svd_u.rightCols(1) = svd_u.rightCols(1) * -1; 
			rotation = svd.matrixV().transpose() * svd_u; // TODO not sure which, U or V or both, need to be transposed
		}
        assert(rotation.determinant() > 0);
		mesh->setRotation(vertexID, rotation);
}

//Calculates the right hand side of the LGS TODO I think maybe there is a bug here
Vector3f calculateB(SimpleMesh *mesh, int i){
    vector<int> neighbors = mesh->getNeighborsOf(i);
    int numNeighbors = neighbors.size();
    Vector3f b = Vector3f::Zero();
    if(mesh->isInFixedVertices(i)){
        if(mesh->isHandle(i)){
            b = mesh->getHandleNewPosition();
        }
        else{
            b = mesh->getVertex(i);
        }
    }
    else{
    
    for ( int neighborID : neighbors)
    {
        b += mesh->getWeight(i, neighborID) * 0.5 * (mesh->getRotation(i)+ mesh->getRotation(neighborID))*(mesh->getVertex(i) - mesh->getVertex(neighborID));
    }
    }
    return b;
}

void estimateVertices(SimpleMesh *mesh){
    MatrixXf b = MatrixXf::Zero(mesh->getNumberOfVertices(), 3);
    for ( int i = 0; i< mesh->getNumberOfVertices(); ++i)
    {
        b.row(i) = calculateB(mesh, i);
    }

    // cout<<"b: \n"<<b<<endl;
    // cout<<"Laplacian: \n"<<mesh->getSystemMatrix()<<endl;
    //Solve LES with Cholesky, L positive definite // TODO test sparse cholesky on sparse eigen matrices
    cout<<"Solving LES ..." <<endl;
    MatrixXf PPrime = mesh->getSystemMatrix().colPivHouseholderQr().solve(b); //Householder should work in any case, later cholesky or something faster
    cout<<"Done!"<<endl;
    cout<<"PPrime Result:" <<endl;
    cout<<PPrime<<endl;
    mesh->setPPrime(PPrime); // set the calculated deformed vertices in the mesh
}

float calculateEnergy(SimpleMesh *mesh){ // TODO not sure if implemented energy function correctly
    float energy= 0.0f;
    for(int i=0; i<mesh->getNumberOfVertices(); ++i){
        vector<int> neighbors = mesh->getNeighborsOf(i);
        // int numNeighbors = neighbors.size();
        float cell_energy=0;

        for ( int j : neighbors)
        {
            Vector3f v= ((mesh->getDeformedVertex(i) - mesh->getDeformedVertex(j)) - mesh->getRotation(i) * (mesh->getVertex(i) - mesh->getVertex(j))); 
            cell_energy += mesh->getWeight(i,j) * v.norm();
        } 
        energy+=cell_energy;
    }
    return energy;
}

void applyDeformation(SimpleMesh *mesh, int handleID, Vector3f handleNewPosition, int iterations){
    mesh->setNewHandlePosition(handleNewPosition);
    // mesh->printPs();
    // mesh->printNewHandlePosition();
    float energy=0.0f;
    int iter=0;
    cout<<"Applying deformation for handle with ID " << handleID << " to new position " << handleNewPosition.x() <<","<< handleNewPosition.y()<< ","<< handleNewPosition.z()<<endl;
    while(iter<iterations){
        cout<<"[Iteration "<<iter<<"]"<<endl;

        // estimateVertices(mesh);
        for(int i=0; i< mesh->getNumberOfVertices(); ++i){
            estimateRotation(mesh, i);
        }
        estimateVertices(mesh);
        float energy_i = calculateEnergy(mesh);        
        cout<< "Iteration: "<< iter<< "  Local error: "<< energy_i << endl;

        mesh->copyPPrime(); // write new vertice locations as current locations in mesh (PPrime -> P)

        energy = energy_i;
        iter++;
    }
    mesh->copyPPrime(); // write new vertice locations as current locations in mesh (PPrime -> P)
    cout << "Resulting energy: "<< energy<< endl; //TODO energy getting bigger rather than smaller :/
    cout << "PPrime[handleID] is "<< mesh->getDeformedVertex(handleID).x() <<","<< mesh->getDeformedVertex(handleID).y()<< ","<< mesh->getDeformedVertex(handleID).z()<<endl;
    // assert(mesh.getDeformedVertex(handleID).x() == handleNewPosition.x());
    // assert(mesh.getDeformedVertex(handleID).y() == handleNewPosition.y());
    // assert(mesh.getDeformedVertex(handleID).z() == handleNewPosition.z());

    // mesh->printPs();
}