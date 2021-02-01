#include "ArapDeformer.h"

#define USE_UNIFORM_WEIGHTS false
#define USE_COTANGENT_WEIGHTS false
#define USE_CONSTANT_WEIGHTS true

#define USE_SPARSE true

ArapDeformer::ArapDeformer() {

}

ArapDeformer::ArapDeformer(SimpleMesh *mesh) {
    m_mesh = *mesh;
    m_num_p = m_mesh.getNumberOfFaces();
    m_num_v = m_mesh.getNumberOfVertices();
    m_cell_rotations = std::vector<MatrixXf>(m_num_v);
    m_b = MatrixXf::Zero(m_num_v, 3);
    m_system_matrix = MatrixXf::Zero(m_num_v, m_num_v);

    for (int i = 0; i < m_num_v; i++) {
        m_cell_rotations[i] = MatrixXf::Zero(3, 3);
    }

    buildWeightMatrix();
}

void ArapDeformer::setHandleConstraint(int handleID, Vector3f newHandlePosition){
    for (int i = 0; i < m_constraints.size(); ++i) {
        if (m_constraints[i].vertexID == handleID) {
            cout<<"Setting handle constraint"<<endl;
            m_constraints[i].position = newHandlePosition;
        }
    }
}


void ArapDeformer::estimateRotation(){
    for(int vertexID=0; vertexID< m_num_v; vertexID++){
        Matrix3f rotation = Matrix3f::Identity();
        vector<int> neighbors = m_mesh.getNeighborsOf(vertexID);
        int numNeighbors = neighbors.size();
        

        MatrixXf PPrime =MatrixXf::Zero(3, numNeighbors);
        MatrixXf P =MatrixXf::Zero(3, numNeighbors);
        MatrixXf D = MatrixXf::Zero(numNeighbors, numNeighbors);

        for ( int j = 0; j< numNeighbors; j++)
        {   
            int neighborVertex = neighbors[j];
            PPrime.col(j) = m_mesh.getDeformedVertex(vertexID) - m_mesh.getDeformedVertex(neighborVertex);
            P.col(j) = m_mesh.getVertex(vertexID) - m_mesh.getVertex(neighborVertex);
            D(j,j) = m_weight_matrix(vertexID, neighborVertex);
        } 

        JacobiSVD<MatrixXf> svd(P * D * PPrime.transpose(), ComputeThinU | ComputeThinV);
        rotation = svd.matrixV() * svd.matrixU().transpose(); 
        if(rotation.determinant() < 0)
        {
            MatrixXf svd_u = svd.matrixU();
            svd_u.rightCols(1) = svd_u.rightCols(1) * -1; 
            rotation = svd.matrixV() * svd_u.transpose(); 
        }
        assert(rotation.determinant() > 0);
        m_cell_rotations[vertexID] =  rotation;
    }
}

Vector3f ArapDeformer::getConstraintI(int id) {
    for (Constraint c : m_constraints) {
        if (c.vertexID == id) {
            return c.position;
        }
    }
    throw std::invalid_argument("received id which is not in fixed vertices");
}

bool ArapDeformer::isInConstraints(int i) {
	for (Constraint c : m_constraints) {
		if (c.vertexID == i) return true;
	}
	return false;
}


void ArapDeformer::updateB(){
    m_b = MatrixXf::Zero(m_num_v, 3);
    for ( int i = 0; i< m_num_v; i++)
    {
        Vector3f sum(0.0f, 0.0f, 0.0f);
        if(isInConstraints(i)){
            sum = getConstraintI(i);
        }
        else{
            vector<int> neighbors = m_mesh.getNeighborsOf(i);
            int numNeighbors = neighbors.size();
        
            for ( int neighborID : neighbors)
            {
                float w_ij = m_weight_matrix(i, neighborID);
                sum += w_ij * 0.5 * (m_cell_rotations[i] + m_cell_rotations[neighborID])*(m_mesh.getVertex(i) - m_mesh.getVertex(neighborID));
            }
        }
        m_b.row(i) = sum;
    }
}



void ArapDeformer::estimateVertices(){
    MatrixXf system_matrix = m_system_matrix;
    cout<<"Solving LES ..." <<endl;
    static JacobiSVD<Eigen::MatrixXf> svd(system_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    MatrixXf result = svd.solve(m_b);
    //cout<<"m_b"<<m_b<<endl;
    cout<<"Done!"<<endl;
    cout<<"Result:" <<endl;
    //cout<<result<<endl;
    m_mesh.setPPrime(result); 
}



float ArapDeformer::calculateEnergy(){
    float energy= 0.0f;
    for(int i=0; i<m_num_v; i++){
        vector<int> neighbors = m_mesh.getNeighborsOf(i);
        float cell_energy=0;

        for ( int j : neighbors)
        {
            Vector3f v= ((m_mesh.getDeformedVertex(i) - m_mesh.getDeformedVertex(j)) - m_cell_rotations[i] * (m_mesh.getVertex(i) - m_mesh.getVertex(j))); 
            cell_energy = m_weight_matrix(i, j) * v.norm();
            energy += cell_energy;
        } 
    }
    return energy;
}

void ArapDeformer::buildWeightMatrix(){
    cout << "Generating Weight Matrix" << endl;
    if(USE_CONSTANT_WEIGHTS){
        cout << "Using constant weights" << endl;
        m_weight_matrix = MatrixXf::Ones(m_num_v, m_num_v);
    }
    else{
        m_weight_matrix = MatrixXf::Zero(m_num_v, m_num_v);
        for (int i = 0; i < m_num_v; i++) {
            if(USE_UNIFORM_WEIGHTS){
                cout << "Using uniform weights" << endl;
                float weight_ij = m_mesh.computeUniformWeightForVertex(i);
                vector<int> neighbors = m_mesh.getNeighborsOf(i);
                for (int j : neighbors){
                    m_weight_matrix(i, j) = weight_ij;
                }
            }
            if (USE_COTANGENT_WEIGHTS){
                cout << "Using cotangent weights" << endl;
                vector<int> neighbors = m_mesh.getNeighborsOf(i);
                for (int j : neighbors){
                    float weight_ij = 0;
                    if (m_weight_matrix(j, i) == 0) 
                        weight_ij = m_mesh.computeCotangentWeightForPair(i, j);
                    else
                        weight_ij = m_weight_matrix(j, i);

                    m_weight_matrix(i, j) = weight_ij;
                }
            }
        }
    }
    
    cout<<"Weight matrix: "<<m_weight_matrix<<endl;
} 


void ArapDeformer::calculateSystemMatrix(){
    for (int i = 0; i < m_num_v; i++) {
        vector<int> neighbors = m_mesh.getNeighborsOf(i);
        int numNeighbors = neighbors.size();
        for (int j = 0; j < numNeighbors; ++j)
        {
            int neighborVertex = neighbors[j];
            m_system_matrix(i, i) += m_weight_matrix(i, neighborVertex);
            m_system_matrix(i, neighborVertex) = -m_weight_matrix(i, neighborVertex);
        }
    }

    for (Constraint c : m_constraints) {
        int i = c.vertexID;
        m_system_matrix.row(i).setZero();
        m_system_matrix(i, i) = 1;
    }
}

void ArapDeformer::initDeformation(vector<int> fixed_points){
    m_constraints.clear();

    for (int i : fixed_points) {
        Constraint c;
        c.vertexID = i;
        c.position = m_mesh.getVertex(i);
        m_constraints.push_back(c);
    }
    calculateSystemMatrix();
}

void ArapDeformer::applyDeformation(int handleID, Vector3f handleNewPosition, int iterations) {

    m_handle_id = handleID;
    m_new_handle_position = handleNewPosition;

    cout << "handleID " << m_handle_id << endl;
    cout << "# fixed Vertices " << m_constraints.size() << endl;

    setHandleConstraint(handleID, handleNewPosition);
    float energy = 999.0f;
    int iter = 0;
    cout << "Applying deformation for handle with ID " << handleID << " to new position " << handleNewPosition.x() << "," << handleNewPosition.y() << "," << handleNewPosition.z() << endl;

    while (iter < iterations && abs(energy) > THRESHOLD) {
        cout << "[Iteration " << iter << "]" << endl;

        estimateRotation();
        updateB();
        estimateVertices();

        float energy_i = calculateEnergy();
        cout << "Iteration: " << iter << "  Local error: " << energy_i << endl;

        m_mesh.copyPPrime();

        energy = energy_i;
        iter++;
    }
    cout << "Resulting energy: " << energy << endl;
    cout << "PPrime[handleID] is " << m_mesh.getDeformedVertex(handleID).x() << "," << m_mesh.getDeformedVertex(handleID).y() << "," << m_mesh.getDeformedVertex(handleID).z() << endl;
    m_mesh.writeMesh("../data/bunny/deformedMesh.off");
}