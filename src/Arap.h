#include <iostream>
#include <fstream>

#include "Eigen.h"
#include "SimpleMesh.h"
using namespace std;

Matrix3f estimateRotation(SimpleMesh mesh, int vertexID) {
		// TODO: Estimate the rotation from source to target points, following the Procrustes algorithm.
		// To compute the singular value decomposition you can use JacobiSVD() from Eigen.
		// Important: The covariance matrices should contain mean-centered source/target points.

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
			Eigen::Matrix3f tmp;
			rotation = svd.matrixV() * Matrix3f::Identity() * Vector3f(1,1,-1) * svd.matrixU().transpose();
		}

		return rotation;
	}

    void applyDeformation(SimpleMesh mesh, int iterations){
        float energy=0.0f;
        cout<<"Applying deformation"<<endl;
        while(iterations>0){
             // initialize b and assign constraints
            //number_of_fixed_verts = len(self.fixed_verts)

            //self.b_array = np.zeros((self.n + number_of_fixed_verts, 3))
            // Constraint b points
            //for i in range(number_of_fixed_verts):
            //  self.b_array[self.n + i] = self.fixed_verts[i][1]

        
            estimateRotation();
            applyRotation(); // estimateVertices()?
            float energy_i = calculateEnergy();
            cout<< "Iterations left: "<< iteration<< "  Local error: "<< energy_i << endl;

            //  if(self.energy_minimized(iteration_energy)):
            //      print("Energy was minimized at iteration", t, " with an energy of ", iteration_energy)
            //      break
            energy = energy_i;
            iteration--;
        }
    }