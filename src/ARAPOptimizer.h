#pragma once

// The Google logging library (GLOG), used in Ceres, has a conflict with Windows defined constants. This definitions prevents GLOG to use the same constants
#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <flann/flann.hpp>

#include "SimpleMesh.h"
#include "NearestNeighbor.h"
#include "PointCloud.h"
#include "ProcrustesAligner.h"
using namespace std;


class ARAPConstraint {
public:
	ARAPConstraint(){}
/*
	template <typename T>
	bool operator()(const T* const pose, T* residual) const {

		// Important: Ceres automatically squares the cost function.
		//loss function set to NULL already takes squared norm of residuals!
		//residual[0] =  0.0f;
		return true;
	}

static ceres::CostFunction* create() {
		return new ceres::AutoDiffCostFunction<ARAPConstraint, 1>(
			new ARAPConstraint()
		);
	}*/

protected:
	const Vector3f m_sourcePoint;
	const Vector3f m_targetPoint;
	const vector<Vertex> m_sourceFan;
	const vector<Vertex> m_targetFan;
};

/**
 * ARAP optimizer - using Ceres for optimization.
 */
class ARAPOptimizer {
public:
	ARAPOptimizer(){}

	Matrix4f deform() {

		// for (int i = 0; i < m_nIterations; ++i) {
		//  	// Prepare constraints
		//  	ceres::Problem problem;
		//  	unsigned nPoints = source.getPoints().size();

		//  	for (unsigned i = 0; i < nPoints; ++i) {
		// 		 const auto& sourcePoint = sourcePoints[i];
		// 		 const auto& targetPoint = targetPoints[i];
		// 		 const auto& sourceFan = sourceFans[i];
		// 		 const auto& targetFan = targetFans[i];

		//  			problem.AddResidualBlock(
		//  				new ceres::AutoDiffCostFunction<ARAPConstraint, 1, 6>(
		//  					new ARAPConstraint(sourcePoint, targetPoint, sourceFan, targetFan)), 
		//  					nullptr,
		//  					poseIncrement.getData());
		//  	}
		// }

		// // Configure options for the solver.
		// ceres::Solver::Options options;
		// configureSolver(options);

		// // Run the solver (for one iteration).
		// ceres::Solver::Summary summary;
		// ceres::Solve(options, &problem, &summary);
		// std::cout << summary.BriefReport() << std::endl;
		// //std::cout << summary.FullReport() << std::endl;

		// // Update the current pose estimate (we always update the pose from the left, using left-increment notation).
		// // Matrix4f matrix = PoseIncrement<double>::convertToMatrix(poseIncrement);
		// // estimatedPose = PoseIncrement<double>::convertToMatrix(poseIncrement) * estimatedPose;
		// // poseIncrement.setZero();

		// std::cout << "Optimization iteration done." << std::endl;
		
		return Matrix4f::Identity();
	}

	void setNumIterations( unsigned n){
		m_nIterations = n;
	}

private:
	unsigned m_nIterations;

	void configureSolver(ceres::Solver::Options& options) {
		// Ceres options.
		options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
		options.use_nonmonotonic_steps = false;
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = 1;
		options.max_num_iterations = 1;
		options.num_threads = 8;
	}


};