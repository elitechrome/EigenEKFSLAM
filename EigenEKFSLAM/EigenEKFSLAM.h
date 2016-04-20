//
//  EigenEKFSLAM.h
//  EigenEKFSLAM
//
//  Created by jaemin on 1/26/16.
//  Copyright © 2016 jaemin. All rights reserved.
//

#ifndef EigenEKFSLAM_h
#define EigenEKFSLAM_h
#include <cmath>

#include "Eigen/Eigen"
#include "Eigen/Sparse"
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Matrix<double, 7, 1> ControlVector;
typedef Eigen::Matrix<double, 8, 1> MeasurementVector;

class EigenEKFSLAM{
private:
    Eigen::VectorXd state;
    Eigen::MatrixXd stateCovariance;
    Eigen::Vector3d stateVehicleT;
    Eigen::Quaterniond stateVehicleR;
    ControlVector controlVector;
    
    //Motion Model, Sensor Model, innovation, kalman gain
    /* Intermediate variables used during filter iteration. */
    Eigen::VectorXd aprioriMean;
    Eigen::MatrixXd aprioriCovariance;
    Eigen::MatrixXd R;
    
public:
    //Mean prediction with deterministic motion model func
    //Covariance prediction with Jacobian of Motion model func
    //Measurement conversion with Sensor model func
    //Calculate the innovation matrix with predicted cov and Jacobian of Sensor model func
    // + Sensor noise model(Constant?)
    //Calculate the Kalman gain with innovation matrix and Jacobian of Sensor model func
    void aprioriEstimate(double dt, ControlVector& c); 
    void measurementEstimate(std::vector<MeasurementVector> &m);
    void addLandmark(MeasurementVector _stateFeature);
    const Eigen::VectorXd getState();
  
public:
    EigenEKFSLAM(Eigen::Vector3d _stateT, Eigen::Quaterniond _stateR);
    ~EigenEKFSLAM(){}
};

#endif /* EigenEKFSLAM_h */
