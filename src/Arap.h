#include <iostream>
#include <fstream>

#include "Eigen.h"
#include "SimpleMesh.h"
using namespace std;

#define THRESHOLD 1.0e-3f

void estimateRotation(SimpleMesh *mesh, int vertexID) {
		// Assume vertices are fix, solve for rotations with Procrustes

        // Sorkine paper: PDP' for vertex i with neighbors j (P contains edges of fan as cols, P' edges of deformed fan, D diagonal with weights w_ij)
		Matrix3f rotation = Matrix3f::Identity();
        vector<int> neighbors = mesh->getNeighborsOf(vertexID);
        int numNeighbors = neighbors.size();

		// // P precomputed in SimpleMesh.loadMesh()
		MatrixXf PPrime =MatrixXf::Zero(3, numNeighbors); // TODO did I miss any initialization of the deformed vertices?
        MatrixXf P =MatrixXf::Zero(3, numNeighbors);
        MatrixXf D = MatrixXf::Zero(numNeighbors, numNeighbors);

		for ( int j = 0; j< numNeighbors; ++j)
		{   
            int neighborVertex = neighbors[j];
			PPrime.col(j) = mesh->getDeformedVertex(vertexID) - mesh->getDeformedVertex(neighborVertex);
            P.col(j) = mesh->getVertex(vertexID) - mesh->getVertex(neighborVertex);
            //D(j,j) = 1.0f;
            D(j,j) = mesh->getWeight(vertexID, neighborVertex);
		} 

        //MatrixXf P = mesh->getPrecomputedP(vertexID);

        //Procrustes
		JacobiSVD<MatrixXf> svd(P * D * PPrime.transpose(), ComputeThinU | ComputeThinV);
		rotation = svd.matrixV() * svd.matrixU().transpose(); // TODO svd() gives A=USV but in paper A=USV' and R=VU' -> both need transpose? Edit: We think that MatrixV from svd function is already V and not V'
		if(rotation.determinant() < 0)
		{
            MatrixXf svd_u = svd.matrixU();
            svd_u.rightCols(1) = svd_u.rightCols(1) * -1; 
			rotation = svd.matrixV() * svd_u.transpose(); // TODO not sure which, U or V or both, need to be transposed Edit: Should be the same calculation as before (V*U')
		}
        assert(rotation.determinant() > 0);
		mesh->setRotation(vertexID, rotation);
}

//Calculates the right hand side of the LGS TODO I think maybe there is a bug here
void calculateB(SimpleMesh *mesh){

    mesh->m_b = MatrixXf::Zero(mesh->getNumberOfVertices(), 3);
    for ( int i = 0; i< mesh->getNumberOfVertices(); ++i)
    {
        vector<int> neighbors = mesh->getNeighborsOf(i);
        int numNeighbors = neighbors.size();
        Vector3f sum(0.0f, 0.0f, 0.0f);
        for ( int neighborID : neighbors)
        {
            float w_ij = mesh->getWeight(i, neighborID); //Added corresponding weights instead of 1.0f
            sum += w_ij * 0.5 * (mesh->getRotation(i)+ mesh->getRotation(neighborID))*(mesh->getVertex(i) - mesh->getVertex(neighborID));
            // b += mesh->getWeight(i, neighborID) * 0.5 * (mesh->getRotation(i)+ mesh->getRotation(neighborID))*(mesh->getVertex(i) - mesh->getVertex(neighborID));
        }

        mesh->m_b.row(i) = sum;
    }
    
}

void estimateVertices(SimpleMesh *mesh){
    // MatrixXf b = MatrixXf::Zero(mesh->getNumberOfVertices(), 3);
    calculateB(mesh);
    MatrixXf systemMatrix = mesh->getSystemMatrix();

    for(int fixedVertex : mesh->getFixedVertices()){

        for(int i=0; i< mesh->getNumberOfVertices(); ++i){
            mesh->m_b.row(i) -= systemMatrix(i, fixedVertex) * mesh->getDeformedVertex(fixedVertex);
			// mesh->m_b.row(i)[1] -= systemMatrix(i, fixedVertex) * mesh->getDeformedVertex(fixedVertex).y();
			// mesh->m_b.row(i)[2] -= systemMatrix(i, fixedVertex) * mesh->getDeformedVertex(fixedVertex).z();
        }

        mesh->m_b.row(fixedVertex) = mesh->getDeformedVertex(fixedVertex);
        for (int i = 0; i < mesh->getNumberOfVertices(); ++i) systemMatrix(fixedVertex, i) = systemMatrix(i, fixedVertex) = 0.0f;
		systemMatrix(fixedVertex, fixedVertex) = 1.0f;
    }
    
    //Solve LES with Cholesky, L positive definite // TODO test sparse cholesky on sparse eigen matrices
    cout<<"Solving LES ..." <<endl;
    // MatrixXf PPrime = mesh->getSystemMatrix().colPivHouseholderQr().solve(b); //Householder should work in any case, later cholesky or something faster
    static JacobiSVD<Eigen::MatrixXf> svd(systemMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    MatrixXf result = svd.solve(mesh->m_b);
    // auto result_x = svd.solve(mesh->m_b.col(0));
    // auto result_y = svd.solve(mesh->m_b.col(1));
    // auto result_z = svd.solve(mesh->m_b.col(2));
    cout<<"m_b"<<mesh->m_b<<endl;

    // MatrixXf result(3, mesh->getNumberOfVertices());
    // result.col(0)=result_x;
    // result.col(1)=result_y;
    // result.col(2)=result_z;
    cout<<"Done!"<<endl;
    cout<<"Result:" <<endl;
    cout<<result<<endl;
    mesh->setPPrime(result); // set the calculated deformed vertices in the mesh
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
            cell_energy += v.norm();
            // cell_energy += mesh->getWeight(i,j) * v.norm();
        } 
        energy+=cell_energy;
    }
    return energy;
}

void applyDeformation(SimpleMesh *mesh, int handleID, Vector3f handleNewPosition, vector<int> fixedPoints, int iterations){
    mesh->applyConstrainedPoints(handleID, handleNewPosition, fixedPoints);
    // mesh->printPs();
    // mesh->printNewHandlePosition();
    float energy=999.0f;
    int iter=0;
    cout<<"Applying deformation for handle with ID " << handleID << " to new position " << handleNewPosition.x() <<","<< handleNewPosition.y()<< ","<< handleNewPosition.z()<<endl;
    while(iter < iterations && energy > THRESHOLD){
        cout<<"[Iteration "<<iter<<"]"<<endl;

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
    cout << "Resulting energy: "<< energy<< endl; //TODO energy getting bigger rather than smaller :/
    cout << "PPrime[handleID] is "<< mesh->getDeformedVertex(handleID).x() <<","<< mesh->getDeformedVertex(handleID).y()<< ","<< mesh->getDeformedVertex(handleID).z()<<endl;
    // assert(mesh.getDeformedVertex(handleID).x() == handleNewPosition.x());
    // assert(mesh.getDeformedVertex(handleID).y() == handleNewPosition.y());
    // assert(mesh.getDeformedVertex(handleID).z() == handleNewPosition.z());

    // mesh->printPs();
}