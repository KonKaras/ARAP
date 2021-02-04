#include "ArapDeformer.h"

typedef Eigen::Triplet<double> T;

ArapDeformer::~ArapDeformer()
{
    m_mesh.writeMesh("../data/bunny/deformdedMesh.off");
}

ArapDeformer::ArapDeformer(){}

ArapDeformer::ArapDeformer(SimpleMesh *mesh, int weight_type, int estimation_type) {
    m_mesh = *mesh;
    m_num_p = m_mesh.getNumberOfFaces();
    m_num_v = m_mesh.getNumberOfVertices();
    m_cell_rotations = std::vector<MatrixXd>(m_num_v);
    m_b = MatrixXd::Zero(m_num_v, 3);    

    for (int i = 0; i < m_num_v; i++) {
        m_cell_rotations[i] = MatrixXd::Zero(3, 3);
    }

    setWeightType(weight_type);
    setDecompositionType(estimation_type);

    buildWeightMatrix();
    calculateSystemMatrix();
}

void ArapDeformer::setWeightType(int weight_type) {
    if (weight_type == 0) {
        use_uniform_weights = true;
        use_constant_weights = false;
        use_cotangent_weights = false;
    }
    else if (weight_type == 1) {
        use_uniform_weights = false;
        use_constant_weights = true;
        use_cotangent_weights = false;
    }
    //if weight_type = 2 or any other number than 1 and 0
    else {
        use_uniform_weights = false;
        use_constant_weights = false;
        use_cotangent_weights = true;
    }
}

void ArapDeformer::setDecompositionType(int estimation_type) {
    /*if (estimation_type == 0) {
        use_simplicial_ldlt = true;
        use_simplicial_llt = false;
        use_sparse_qr = false;
        use_sparse_lu = false;
    }
    else if (estimation_type == 1) {
        use_simplicial_ldlt = false;
        use_simplicial_llt = true;
        use_sparse_qr = false;
        use_sparse_lu = false;
    }
    else */
    if (estimation_type == 0) {
        use_simplicial_ldlt = false;
        use_simplicial_llt = false;
        use_sparse_qr = true;
        use_sparse_lu = false;
    }
    else if(estimation_type == 1) {
        use_simplicial_ldlt = false;
        use_simplicial_llt = false;
        use_sparse_qr = false;
        use_sparse_lu = true;
    }
    else if(estimation_type == 2){
        use_simplicial_ldlt = false;
        use_simplicial_llt = false;
        use_sparse_qr = false;
        use_sparse_lu = false;
    }

    if (use_simplicial_ldlt || use_simplicial_llt || use_sparse_qr || use_sparse_lu) {
        use_sparse_matrices = true;
    }
    else {
        use_sparse_matrices = false;
    }
}

void ArapDeformer::setHandleConstraint(int handleID, Vector3d newHandlePosition){
    for (int i = 0; i < m_constraints.size(); i++) {
        if (m_constraints[i].vertexID == handleID) {
            //std::cout<<"Setting handle constraint"<<endl;
            m_constraints[i].position = newHandlePosition;
        }
    }
}


void ArapDeformer::estimateRotation(){
    for(int vertexID=0; vertexID< m_num_v; vertexID++){
        Matrix3d rotation = Matrix3d::Identity();
        vector<int> neighbors = m_mesh.getNeighborsOf(vertexID);
        int numNeighbors = neighbors.size();
        

        MatrixXd PPrime =MatrixXd::Zero(3, numNeighbors);
        MatrixXd P =MatrixXd::Zero(3, numNeighbors);
        MatrixXd D = MatrixXd::Zero(numNeighbors, numNeighbors);

        for ( int j = 0; j< numNeighbors; j++)
        {   
            int neighborVertex = neighbors[j];
            PPrime.col(j) = m_mesh.getDeformedVertex(vertexID) - m_mesh.getDeformedVertex(neighborVertex);
            P.col(j) = m_mesh.GetVertexOriginal(vertexID) - m_mesh.GetVertexOriginal(neighborVertex);
            D(j,j) = m_weight_matrix(vertexID, neighborVertex);
        } 

        JacobiSVD<MatrixXd> svd(P * D * PPrime.transpose(), ComputeThinU | ComputeThinV);
        rotation = svd.matrixV() * svd.matrixU().transpose(); 
        if(rotation.determinant() < 0)
        {
            MatrixXd svd_u = svd.matrixU();
            svd_u.rightCols(1) = svd_u.rightCols(1) * -1; 
            rotation = svd.matrixV() * svd_u.transpose(); 
        }
        assert(rotation.determinant() > 0);
        m_cell_rotations[vertexID] =  rotation;
    }
}

Vector3d ArapDeformer::getConstraintI(int id) {
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
    m_b = MatrixXd::Zero(m_num_v, 3);
    for ( int i = 0; i< m_num_v; i++)
    {
        Vector3d sum(0.0f, 0.0f, 0.0f);
        if(isInConstraints(i)){
            sum = getConstraintI(i);
        }
        else{
            vector<int> neighbors = m_mesh.getNeighborsOf(i);
            int numNeighbors = neighbors.size();
        
            for ( int neighborID : neighbors)
            {
                double w_ij = m_weight_matrix(i, neighborID);
                sum += w_ij * 0.5 * (m_cell_rotations[i] + m_cell_rotations[neighborID])*(m_mesh.GetVertexOriginal(i) - m_mesh.GetVertexOriginal(neighborID));
            }
        }
        m_b.row(i) = sum;
    }
}



void ArapDeformer::estimateVertices(){
    std::cout<<"Solving LES ..." <<endl;
    MatrixXd result;
    if(!use_sparse_matrices){
         MatrixXd system_matrix = m_system_matrix;
         static JacobiSVD<Eigen::MatrixXd> svd(system_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
         result = svd.solve(m_b);
         m_mesh.setPPrime(result);
    }
    else{
        /*if(use_simplicial_ldlt){
            SimplicialLDLT<SparseMatrix<double>> solver;
            solver.compute(m_system_matrix_sparse);
            if(solver.info()!=Success) {
                cout<<"Decomposition failed!"<<endl;
                return;
            }
            result = solver.solve(m_b);
            if(solver.info()!=Success) {
                cout<<"Solving failed"<<endl;
                return;
            }
            //std::cout<<"Result:" <<endl;
            //std::cout<<result<<endl;
            m_mesh.setPPrime(result); 
        }
        else if(use_simplicial_llt){
            SimplicialLLT<SparseMatrix<double>, Lower, NaturalOrdering<int>> solver;
            solver.compute(m_system_matrix_sparse);
            if(solver.info()!=Success) {
                cout<<"Decomposition failed!"<<endl;
                return;
            }
            result = solver.solve(m_b);
            if(solver.info()!=Success) {
                cout<<"Solving failed"<<endl;
                return;
            }
            m_mesh.setPPrime(result);
        }
        else */
        if(use_sparse_lu){
            SparseLU<SparseMatrix<double>> solver;
            solver.compute(m_system_matrix_sparse);
            if(solver.info()!=Success) {
                cout<<"Decomposition failed!"<<endl;
                return;
            }
            result = solver.solve(m_b);
            if(solver.info()!=Success) {
                cout<<"Solving failed"<<endl;
                return;
            }
            m_mesh.setPPrime(result);
        }
        else if (use_sparse_qr){
            SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
            m_system_matrix_sparse.makeCompressed();
            solver.compute(m_system_matrix_sparse);
            if(solver.info()!=Success) {
                cout<<"Decomposition failed!"<<endl;
                return;
            }
            result = solver.solve(m_b);
            if(solver.info()!=Success) {
                cout<<"Solving failed"<<endl;
                return;
            }
            m_mesh.setPPrime(result);
        }
        else{
            cout<<"No decomposition chosen!"<<endl;
            return;
        }
    }
    std::cout<<"Done!"<<endl;
}



double ArapDeformer::calculateEnergy(){
    double energy= 0.0f;
    for(int i=0; i<m_num_v; i++){
        vector<int> neighbors = m_mesh.getNeighborsOf(i);
        double cell_energy=0;

        for ( int j : neighbors)
        {
            Vector3d v = ((m_mesh.getDeformedVertex(i) - m_mesh.getDeformedVertex(j)) - m_cell_rotations[i] * (m_mesh.GetVertexOriginal(i) - m_mesh.GetVertexOriginal(j)));//(m_mesh.getVertex(i) - m_mesh.getVertex(j))); 
            cell_energy += m_weight_matrix(i, j) * v.squaredNorm();
        } 
        energy += m_weight_matrix(i, i) * cell_energy;   
    }
    return energy;
}

void ArapDeformer::buildWeightMatrix(){
    std::cout << "Generating Weight Matrix" << endl;
    if(use_constant_weights){
        std::cout << "Using constant weights" << endl;
        m_weight_matrix = MatrixXd::Ones(m_num_v, m_num_v);
    }
    else{
        m_weight_matrix = MatrixXd::Zero(m_num_v, m_num_v);
        for (int i = 0; i < m_num_v; i++) {
            if(use_uniform_weights){
                std::cout << "Using uniform weights" << endl;
                double weight_ij = m_mesh.computeUniformWeightForVertex(i);
                vector<int> neighbors = m_mesh.getNeighborsOf(i);
                for (int j : neighbors){
                    m_weight_matrix(i, j) = weight_ij;
                }
            }
            if (use_cotangent_weights){
                std::cout << "Using cotangent weights" << endl;
                vector<int> neighbors = m_mesh.getNeighborsOf(i);
                for (int j : neighbors){
                    double weight_ij = 0;
                    if (m_weight_matrix(j, i) == 0) 
                        weight_ij = m_mesh.computeCotangentWeightForPair(i, j);
                    else
                        weight_ij = m_weight_matrix(j, i);

                    m_weight_matrix(i, j) = weight_ij;
                }
            }
            m_weight_matrix(i, i) = 1;
        }
    }
    
    //std::cout<<"Weight matrix: "<<m_weight_matrix<<endl;
} 


void ArapDeformer::calculateSystemMatrix(){
    m_system_matrix = MatrixXd::Zero(m_num_v, m_num_v);
    for (int i = 0; i < m_num_v; i++) {
        vector<int> neighbors = m_mesh.getNeighborsOf(i);
        int numNeighbors = neighbors.size();
        for (int j = 0; j < numNeighbors; ++j){
            int neighborVertex = neighbors[j];
            m_system_matrix(i, i) += m_weight_matrix(i, neighborVertex);
            m_system_matrix(i, neighborVertex) = -m_weight_matrix(i, neighborVertex);
        }
    }
   
    // for (Constraint c : m_constraints) {
    //     int i = c.vertexID;
    //     m_system_matrix.row(i).setZero();
    //     m_system_matrix(i, i) = 1;
    // }

    m_system_matrix_original = m_system_matrix;

    // if (USE_SPARSE_MATRICES)
    //     m_system_matrix_sparse = m_system_matrix.sparseView();

}

void ArapDeformer::updateSystemMatrix(){

	m_system_matrix = m_system_matrix_original;
	for (Constraint c : m_constraints) {
		int i = c.vertexID;
		m_system_matrix.row(i).setZero();
		//m_system_matrix.col(i).setZero();
		m_system_matrix(i, i) = 1;
	}

    if (use_sparse_matrices)
        m_system_matrix_sparse = m_system_matrix.sparseView();
}

void ArapDeformer::initDeformation(vector<int> fixed_points){
    m_constraints.clear();
    for (int i : fixed_points) {
        Constraint c;
        c.vertexID = i;
        c.position = m_mesh.GetVertexOriginal(i);
        m_constraints.push_back(c);
    }
    updateSystemMatrix();

    calculateSystemMatrix();
}

void ArapDeformer::applyDeformation(vector<int> fixed_points, int handleID, Vector3d handleNewPosition, int iterations) {
    m_constraints.clear();

    for (int i : fixed_points) {
        Constraint c;
        c.vertexID = i;
        c.position = m_mesh.GetVertexOriginal(i);
        m_constraints.push_back(c);
    }

    updateSystemMatrix();

    m_handle_id = handleID;
    m_new_handle_position = handleNewPosition;
    // calculateSystemMatrix();

    std::cout << "handleID " << m_handle_id << endl;
    std::cout << "# fixed Vertices " << m_constraints.size() << endl;
    std::cout <<"Using sparse matrices: "<<use_sparse_matrices<<endl;
    //std::cout << "NonZeros in sparse matrix: "<<m_system_matrix_sparse.nonZeros()<<endl;

    setHandleConstraint(handleID, handleNewPosition);
    double energy = 0;
    double energy_prev = 1;
    int iter = 0;
    cout << "Applying deformation for handle with ID " << handleID << " at position " << m_mesh.getVertex(handleID).transpose()<< " to new position " << handleNewPosition.x() << "," << handleNewPosition.y() << "," << handleNewPosition.z() << endl;

    while (iter < iterations) {
        cout << "[Iteration " << iter << "]" << endl;

        estimateRotation();
        updateB();
        estimateVertices();

        double energy_i = calculateEnergy();        
        std::cout<< "Iteration: "<< iter<< "  Local error: "<< energy_i << endl;

        m_mesh.copyPPrime();

        energy = energy_i;

        if (abs(energy_prev - energy) < THRESHOLD) break;

        energy_prev = energy;
        iter++;
    }
    std::cout << "Resulting energy: "<< energy<< endl; 
    std::cout << "PPrime[handleID] is "<< m_mesh.getDeformedVertex(handleID).x() <<","<< m_mesh.getDeformedVertex(handleID).y()<< ","<< m_mesh.getDeformedVertex(handleID).z()<<endl;
    //m_mesh.writeMesh("../data/bunny/deformedMesh.off"); 
}