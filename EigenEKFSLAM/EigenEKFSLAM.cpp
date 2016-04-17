//
//  EigenEKFSLAM.cpp
//  EigenEKFSLAM
//
//  Created by jaemin on 1/26/16.
//  Copyright © 2016 jaemin. All rights reserved.
//
#include "EigenEKFSLAM.h"
EigenEKFSLAM::EigenEKFSLAM(Eigen::Vector3d _stateT, Eigen::Quaterniond _stateR):state(7), stateCovariance(7,7)
{
    
    Eigen::VectorXd state;
    Eigen::MatrixXd stateCovariance;
    
    Eigen::Vector3d stateVehicleT;
    Eigen::Quaterniond stateVehicleR;
    ControlVector controlVector;
    
    //Motion Model, Sensor Model, innovation, kalman gain
    /* Intermediate variables used during filter iteration. */
    
    Eigen::VectorXd apriori_mean;
    Eigen::MatrixXd apriori_covariance;
}
void EigenEKFSLAM::aprioriEstimate(double dt, ControlVector& c)
{
    //ToDo: modify code with
    double qw = state[3], qx = state[4], qy = state[5], qz = state[6];
    state[0] += c[0]*(qw*qw +qx*qx -qy*qy -qz*qz)   + 2*c[1]*(qx*qy -qw*qz)         + 2*c[2]*(qx*qz +qw*qy);
    state[1] += 2*c[0]*(qx*qy +qw*qz)               + (qw*qw -qx*qx +qy*qy -qz*qz)  + 2*c[2]*(qy*qz -qw*qx);
    state[2] += 2*c[0]*(qx*qz -qw*qy)               + 2*c[1]*(qy*qz +qw*qx)         + c[2]*(qw*qw -qx*qx -qy*qy +qz*qz);
    state[3] = qw*c[3]-qx*c[5]-qy*c[4]-qz*c[6];
    state[4] = qx*c[3]+qw*c[5]-qz*c[4]+qy*c[6];
    state[5] = qy*c[3]+qz*c[5]+qw*c[4]-qx*c[6];
    state[6] = qz*c[3]-qy*c[5]+qx*c[4]+qw*c[6];
    
    Eigen::MatrixXd vehicleCov = stateCovariance.block(0, 0, 7, 7);
    Eigen::MatrixXd noiseCov(7, 7);
    Eigen::MatrixXd J_fv_xv(7, 7);
    Eigen::MatrixXd J_u_xv(7, 7);
    
    //ToDo: Determine control noise cofficients
    noiseCov <<
    0.5*c[0], 0, 0, 0, 0, 0, 0,
    0, 0.5*c[1], 0, 0, 0, 0, 0,
    0, 0, 0.5*c[2], 0, 0, 0, 0,
    0, 0, 0, 0.5*c[3], 0, 0, 0,
    0, 0, 0, 0, 0.5*c[4], 0, 0,
    0, 0, 0, 0, 0, 0.5*c[5], 0,
    0, 0, 0, 0, 0, 0, 0.5*c[6];
    
    //Jacobian of predicted vehicle state and vehicle state
    Eigen::MatrixXd A(3,4), B(4,4), C(3,3), D(4,4);
    A <<
    qw*c[0]-qz*c[1]+qy*c[2], qx*c[0]+qy*c[1]+qz*c[2], -qy*c[0]+qx*c[1]+qw*c[2], -qz*c[0]-qw*c[1]+qz*c[2],
    qz*c[0]+qw*c[1]-qx*c[2], qy*c[0]-qx*c[1]-qw*c[2], qx*c[0]+qy*c[1]+qz*c[2], qx*c[0]-qz*c[1]+qy*c[2],
    -qy*c[0]+qx*c[1]+qw*c[2], qz*c[0]+qw*c[1]-qx*c[2], -qw*c[0]+qz*c[1]-qy*c[2], qx*c[0]-qy*c[1]-qz*c[2];
    
    A *= 2;
    
    B <<
    c[3], -c[4], -c[5], -c[6],
    c[4], c[3], c[6], -c[5],
    c[5], -c[6], c[3], c[4],
    c[6], c[5], -c[4], c[3];
    
    C <<
    qw*qw +qx*qx -qy*qy -qz*qz, 2*(qx*qy-qw*qz), 2*(qw*qy +qx*qz),
    2*(qw*qz+qx*qy), qw*qw -qx*qx +qy*qy -qz*qz, 2*(qy*qz -qw*qx),
    2*(qx*qz-qw*qy), 2*(qw*qx+qy*qz), qw*qw -qx*qx -qy*qy +qz*qz;
    
    D <<
    qx, -qx, -qy, -qz,
    qx, qw, -qz, qy,
    qy, qz,  qw, -qx,
    qz, -qz, qy, qw;
    
    
    J_fv_xv.block(0, 0, 3, 3).setIdentity();
    J_fv_xv.block(3, 0, 4, 3).setZero();
    J_fv_xv.block(0, 3, 3, 4) = A;
    J_fv_xv.block(3, 3, 4, 4) = B;
    
    J_u_xv.block(0, 0, 3, 3) = C;
    J_u_xv.block(0, 3, 3, 3).setZero();
    J_u_xv.block(3, 0, 4, 3).setZero();
    J_u_xv.block(3, 3, 4, 4) = D;
    
    stateCovariance.block(0, 0, 7, 7) = J_fv_xv * vehicleCov * J_fv_xv.transpose() + J_u_xv * noiseCov * J_u_xv.transpose();
    
}
void EigenEKFSLAM::measurementEstimate(std::vector<MeasurementVector> &m)
{
    Eigen::MatrixXd J_h_x(stateCovariance.rows()-state.rows(),stateCovariance.cols());
    J_h_x.setZero();
    Eigen::VectorXd y_k(state.size());
    Eigen::MatrixXd J_hi_xv(7, 7);
    Eigen::MatrixXd J_hi_yi(7, 7);
    double x0 = state[0], y0 = state[1], z0 = state[2], qw0 = state[3], qx0 = state[4], qy0 = state[5], qz0 = state[6];

    //for each observed landmark
    for (int i = 0; i < m.size(); i++){
        int s = lround(m[i][7]);
        if (s > m.size()) {
            addLandmark(m[i]);
            continue;
        }
        //landmakr absolute 3D point
        double xi = m[s][0], yi = m[s][1], zi = m[s][2], qwi = m[s][3], qxi = m[s][4], qyi = m[s][5], qzi = m[s][6];
    
    J_hi_xv
    <<
        - (2*qw0*qy0 + 2*qx0*qz0)*(2*qwi*qyi + 2*qxi*qzi) - (2*qw0*qz0 - 2*qx0*qy0)*(2*qwi*qzi - 2*qxi*qyi) - (2*pow(qy0,2) + 2*pow(qz0,2) - 1)*(2*pow(qyi,2) + 2*pow(qzi,2) - 1),   (2*qw0*qz0 + 2*qx0*qy0)*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) - (2*qwi*qzi - 2*qxi*qyi)*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + (2*qw0*qx0 - 2*qy0*qz0)*(2*qwi*qyi + 2*qxi*qzi),   (2*qwi*qyi + 2*qxi*qzi)*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) - (2*qw0*qy0 - 2*qx0*qz0)*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + (2*qw0*qx0 + 2*qy0*qz0)*(2*qwi*qzi - 2*qxi*qyi), y0*(2*qz0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 2*qx0*(2*qwi*qyi + 2*qxi*qzi)) - z0*(2*qy0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) - 2*qx0*(2*qwi*qzi - 2*qxi*qyi)) - x0*(2*qy0*(2*qwi*qyi + 2*qxi*qzi) + 2*qz0*(2*qwi*qzi - 2*qxi*qyi)), y0*(2*qy0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qyi + 2*qxi*qzi) - 4*qx0*(2*qwi*qzi - 2*qxi*qyi)) + z0*(2*qz0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qzi - 2*qxi*qyi) + 4*qx0*(2*qwi*qyi + 2*qxi*qzi)) + x0*(2*qy0*(2*qwi*qzi - 2*qxi*qyi) - 2*qz0*(2*qwi*qyi + 2*qxi*qzi)), z0*(4*qy0*(2*qwi*qyi + 2*qxi*qzi) - 2*qw0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 2*qz0*(2*qwi*qzi - 2*qxi*qyi)) - x0*(4*qy0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qyi + 2*qxi*qzi) - 2*qx0*(2*qwi*qzi - 2*qxi*qyi)) + y0*(2*qx0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) - 2*qz0*(2*qwi*qyi + 2*qxi*qzi)), z0*(2*qx0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 2*qy0*(2*qwi*qzi - 2*qxi*qyi)) - y0*(2*qy0*(2*qwi*qyi + 2*qxi*qzi) - 2*qw0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 4*qz0*(2*qwi*qzi - 2*qxi*qyi)) - x0*(4*qz0*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qzi - 2*qxi*qyi) + 2*qx0*(2*qwi*qyi + 2*qxi*qzi)),
        (2*qwi*qzi + 2*qxi*qyi)*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) - (2*qw0*qz0 - 2*qx0*qy0)*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + (2*qw0*qy0 + 2*qx0*qz0)*(2*qwi*qxi - 2*qyi*qzi), - (2*qw0*qx0 - 2*qy0*qz0)*(2*qwi*qxi - 2*qyi*qzi) - (2*qw0*qz0 + 2*qx0*qy0)*(2*qwi*qzi + 2*qxi*qyi) - (2*pow(qx0,2) + 2*pow(qz0,2) - 1)*(2*pow(qxi,2) + 2*pow(qzi,2) - 1),   (2*qw0*qx0 + 2*qy0*qz0)*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) - (2*qwi*qxi - 2*qyi*qzi)*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + (2*qw0*qy0 - 2*qx0*qz0)*(2*qwi*qzi + 2*qxi*qyi), z0*(2*qx0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 2*qy0*(2*qwi*qzi + 2*qxi*qyi)) - x0*(2*qz0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) - 2*qy0*(2*qwi*qxi - 2*qyi*qzi)) - y0*(2*qx0*(2*qwi*qxi - 2*qyi*qzi) + 2*qz0*(2*qwi*qzi + 2*qxi*qyi)), x0*(2*qy0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 2*qz0*(2*qwi*qxi - 2*qyi*qzi)) - z0*(4*qx0*(2*qwi*qxi - 2*qyi*qzi) - 2*qw0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 2*qz0*(2*qwi*qzi + 2*qxi*qyi)) - y0*(4*qx0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qxi - 2*qyi*qzi) + 2*qy0*(2*qwi*qzi + 2*qxi*qyi)), x0*(2*qx0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qxi - 2*qyi*qzi) + 4*qy0*(2*qwi*qzi + 2*qxi*qyi)) + z0*(2*qz0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qzi + 2*qxi*qyi) - 4*qy0*(2*qwi*qxi - 2*qyi*qzi)) - y0*(2*qx0*(2*qwi*qzi + 2*qxi*qyi) - 2*qz0*(2*qwi*qxi - 2*qyi*qzi)), x0*(2*qx0*(2*qwi*qxi - 2*qyi*qzi) - 2*qw0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 4*qz0*(2*qwi*qzi + 2*qxi*qyi)) - y0*(4*qz0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + 2*qw0*(2*qwi*qzi + 2*qxi*qyi) - 2*qy0*(2*qwi*qxi - 2*qyi*qzi)) + z0*(2*qy0*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) - 2*qx0*(2*qwi*qzi + 2*qxi*qyi)),
        (2*qw0*qy0 + 2*qx0*qz0)*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) - (2*qwi*qyi - 2*qxi*qzi)*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + (2*qw0*qz0 - 2*qx0*qy0)*(2*qwi*qxi + 2*qyi*qzi),   (2*qwi*qxi + 2*qyi*qzi)*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) - (2*qw0*qx0 - 2*qy0*qz0)*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + (2*qw0*qz0 + 2*qx0*qy0)*(2*qwi*qyi - 2*qxi*qzi), - (2*qw0*qx0 + 2*qy0*qz0)*(2*qwi*qxi + 2*qyi*qzi) - (2*qw0*qy0 - 2*qx0*qz0)*(2*qwi*qyi - 2*qxi*qzi) - (2*pow(qx0,2) + 2*pow(qy0,2) - 1)*(2*pow(qxi,2) + 2*pow(qyi,2) - 1), x0*(2*qy0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 2*qz0*(2*qwi*qxi + 2*qyi*qzi)) - y0*(2*qx0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) - 2*qz0*(2*qwi*qyi - 2*qxi*qzi)) - z0*(2*qx0*(2*qwi*qxi + 2*qyi*qzi) + 2*qy0*(2*qwi*qyi - 2*qxi*qzi)), y0*(4*qx0*(2*qwi*qxi + 2*qyi*qzi) - 2*qw0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 2*qy0*(2*qwi*qyi - 2*qxi*qzi)) - z0*(4*qx0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 2*qw0*(2*qwi*qxi + 2*qyi*qzi) - 2*qz0*(2*qwi*qyi - 2*qxi*qzi)) + x0*(2*qz0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) - 2*qy0*(2*qwi*qxi + 2*qyi*qzi)), y0*(2*qz0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 2*qx0*(2*qwi*qyi - 2*qxi*qzi)) - z0*(4*qy0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 2*qw0*(2*qwi*qyi - 2*qxi*qzi) + 2*qz0*(2*qwi*qxi + 2*qyi*qzi)) - x0*(2*qx0*(2*qwi*qxi + 2*qyi*qzi) - 2*qw0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 4*qy0*(2*qwi*qyi - 2*qxi*qzi)), x0*(2*qx0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 2*qw0*(2*qwi*qxi + 2*qyi*qzi) - 4*qz0*(2*qwi*qyi - 2*qxi*qzi)) + y0*(2*qy0*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + 2*qw0*(2*qwi*qyi - 2*qxi*qzi) + 4*qz0*(2*qwi*qxi + 2*qyi*qzi)) + z0*(2*qx0*(2*qwi*qyi - 2*qxi*qzi) - 2*qy0*(2*qwi*qxi + 2*qyi*qzi)),
        0,0,0,qwi,qxi,qyi,qzi,
        0,0,0,qxi,-qwi,-qzi,qyi,
        0,0,0,qyi,qzi,-qwi,-qxi,
        0,0,0,qzi,-qyi,qxi,-qwi;
        

    
    J_hi_yi <<
        1, 0, 0, z0*(2*qyi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 2*qzi*(2*qw0*qx0 + 2*qy0*qz0)) - y0*(2*qzi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) - 2*qyi*(2*qw0*qx0 - 2*qy0*qz0)) - x0*(2*qyi*(2*qw0*qy0 + 2*qx0*qz0) + 2*qzi*(2*qw0*qz0 - 2*qx0*qy0)),                                                                                                 y0*(2*qyi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 2*qzi*(2*qw0*qx0 - 2*qy0*qz0)) + z0*(2*qzi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) - 2*qyi*(2*qw0*qx0 + 2*qy0*qz0)) + x0*(2*qyi*(2*qw0*qz0 - 2*qx0*qy0) - 2*qzi*(2*qw0*qy0 + 2*qx0*qz0)), y0*(2*qxi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qx0 - 2*qy0*qz0) + 4*qyi*(2*qw0*qz0 + 2*qx0*qy0)) - x0*(4*qyi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qy0 + 2*qx0*qz0) - 2*qxi*(2*qw0*qz0 - 2*qx0*qy0)) - z0*(2*qxi*(2*qw0*qx0 + 2*qy0*qz0) - 2*qwi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 4*qyi*(2*qw0*qy0 - 2*qx0*qz0)), y0*(2*qxi*(2*qw0*qx0 - 2*qy0*qz0) - 2*qwi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 4*qzi*(2*qw0*qz0 + 2*qx0*qy0)) - x0*(4*qzi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qz0 - 2*qx0*qy0) + 2*qxi*(2*qw0*qy0 + 2*qx0*qz0)) + z0*(2*qxi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 2*qwi*(2*qw0*qx0 + 2*qy0*qz0) - 4*qzi*(2*qw0*qy0 - 2*qx0*qz0)),
        0, 1, 0, x0*(2*qzi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 2*qxi*(2*qw0*qy0 + 2*qx0*qz0)) - z0*(2*qxi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) - 2*qzi*(2*qw0*qy0 - 2*qx0*qz0)) - y0*(2*qxi*(2*qw0*qx0 - 2*qy0*qz0) + 2*qzi*(2*qw0*qz0 + 2*qx0*qy0)), x0*(2*qyi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qy0 + 2*qx0*qz0) - 4*qxi*(2*qw0*qz0 - 2*qx0*qy0)) - y0*(4*qxi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qx0 - 2*qy0*qz0) + 2*qyi*(2*qw0*qz0 + 2*qx0*qy0)) + z0*(4*qxi*(2*qw0*qx0 + 2*qy0*qz0) - 2*qwi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 2*qyi*(2*qw0*qy0 - 2*qx0*qz0)),                                                                                                 x0*(2*qxi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) - 2*qzi*(2*qw0*qy0 + 2*qx0*qz0)) + z0*(2*qzi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 2*qxi*(2*qw0*qy0 - 2*qx0*qz0)) - y0*(2*qxi*(2*qw0*qz0 + 2*qx0*qy0) - 2*qzi*(2*qw0*qx0 - 2*qy0*qz0)), z0*(2*qyi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 2*qwi*(2*qw0*qy0 - 2*qx0*qz0) + 4*qzi*(2*qw0*qx0 + 2*qy0*qz0)) - y0*(4*qzi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qz0 + 2*qx0*qy0) - 2*qyi*(2*qw0*qx0 - 2*qy0*qz0)) - x0*(2*qyi*(2*qw0*qy0 + 2*qx0*qz0) - 2*qwi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 4*qzi*(2*qw0*qz0 - 2*qx0*qy0)),
        0, 0, 1, y0*(2*qxi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 2*qyi*(2*qw0*qz0 + 2*qx0*qy0)) - x0*(2*qyi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) - 2*qxi*(2*qw0*qz0 - 2*qx0*qy0)) - z0*(2*qxi*(2*qw0*qx0 + 2*qy0*qz0) + 2*qyi*(2*qw0*qy0 - 2*qx0*qz0)), x0*(2*qzi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qz0 - 2*qx0*qy0) + 4*qxi*(2*qw0*qy0 + 2*qx0*qz0)) - y0*(4*qxi*(2*qw0*qx0 - 2*qy0*qz0) - 2*qwi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 2*qzi*(2*qw0*qz0 + 2*qx0*qy0)) - z0*(4*qxi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 2*qwi*(2*qw0*qx0 + 2*qy0*qz0) - 2*qzi*(2*qw0*qy0 - 2*qx0*qz0)), x0*(4*qyi*(2*qw0*qy0 + 2*qx0*qz0) - 2*qwi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 2*qzi*(2*qw0*qz0 - 2*qx0*qy0)) + y0*(2*qzi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + 2*qwi*(2*qw0*qz0 + 2*qx0*qy0) - 4*qyi*(2*qw0*qx0 - 2*qy0*qz0)) - z0*(4*qyi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + 2*qwi*(2*qw0*qy0 - 2*qx0*qz0) + 2*qzi*(2*qw0*qx0 + 2*qy0*qz0)),                                                                                                 x0*(2*qxi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + 2*qyi*(2*qw0*qz0 - 2*qx0*qy0)) + y0*(2*qyi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) - 2*qxi*(2*qw0*qz0 + 2*qx0*qy0)) + z0*(2*qxi*(2*qw0*qy0 - 2*qx0*qz0) - 2*qyi*(2*qw0*qx0 + 2*qy0*qz0)),
        0, 0, 0,qw0,qx0,qy0,qz0,
        0, 0, 0,-qx0,qw0,qz0,-qy0,
        0, 0, 0,-qy0,-qz0,qw0,qx0,
        0, 0, 0,-qz0,qy0,-qx0,qw0;
        
        Eigen::VectorXd h;
        h << xi - x0*((2*qw0*qy0 + 2*qx0*qz0)*(2*qwi*qyi + 2*qxi*qzi) + (2*qw0*qz0 - 2*qx0*qy0)*(2*qwi*qzi - 2*qxi*qyi) + (2*pow(qy0,2) + 2*pow(qz0,2) - 1)*(2*pow(qyi,2) + 2*pow(qzi,2) - 1)) + y0*((2*qw0*qz0 + 2*qx0*qy0)*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) - (2*qwi*qzi - 2*qxi*qyi)*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + (2*qw0*qx0 - 2*qy0*qz0)*(2*qwi*qyi + 2*qxi*qzi)) + z0*((2*qwi*qyi + 2*qxi*qzi)*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) - (2*qw0*qy0 - 2*qx0*qz0)*(2*pow(qyi,2) + 2*pow(qzi,2) - 1) + (2*qw0*qx0 + 2*qy0*qz0)*(2*qwi*qzi - 2*qxi*qyi)), yi + x0*((2*qwi*qzi + 2*qxi*qyi)*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) - (2*qw0*qz0 - 2*qx0*qy0)*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) + (2*qw0*qy0 + 2*qx0*qz0)*(2*qwi*qxi - 2*qyi*qzi)) - y0*((2*qw0*qx0 - 2*qy0*qz0)*(2*qwi*qxi - 2*qyi*qzi) + (2*qw0*qz0 + 2*qx0*qy0)*(2*qwi*qzi + 2*qxi*qyi) + (2*pow(qx0,2) + 2*pow(qz0,2) - 1)*(2*pow(qxi,2) + 2*pow(qzi,2) - 1)) + z0*((2*qw0*qx0 + 2*qy0*qz0)*(2*pow(qxi,2) + 2*pow(qzi,2) - 1) - (2*qwi*qxi - 2*qyi*qzi)*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) + (2*qw0*qy0 - 2*qx0*qz0)*(2*qwi*qzi + 2*qxi*qyi)), zi + x0*((2*qw0*qy0 + 2*qx0*qz0)*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) - (2*qwi*qyi - 2*qxi*qzi)*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) + (2*qw0*qz0 - 2*qx0*qy0)*(2*qwi*qxi + 2*qyi*qzi)) + y0*((2*qwi*qxi + 2*qyi*qzi)*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) - (2*qw0*qx0 - 2*qy0*qz0)*(2*pow(qxi,2) + 2*pow(qyi,2) - 1) + (2*qw0*qz0 + 2*qx0*qy0)*(2*qwi*qyi - 2*qxi*qzi)) - z0*((2*qw0*qx0 + 2*qy0*qz0)*(2*qwi*qxi + 2*qyi*qzi) + (2*qw0*qy0 - 2*qx0*qz0)*(2*qwi*qyi - 2*qxi*qzi) + (2*pow(qx0,2) + 2*pow(qy0,2) - 1)*(2*pow(qxi,2) + 2*pow(qyi,2) - 1)), qw0*qwi + qx0*qxi + qy0*qyi + qz0*qzi, qw0*qxi - qwi*qx0 - qy0*qzi + qyi*qz0, qw0*qyi - qwi*qy0 + qx0*qzi - qxi*qz0, qw0*qzi - qwi*qz0 - qx0*qyi + qxi*qy0;

        y_k.block(7*s, 0, 7, 1) = m[s].block(0, 0, 7, 1) - h;
        J_h_x.block(0, 7*s, 7, 7) = J_hi_xv;
        J_h_x.block(7*(s+1), 7*s, 7, 7) = J_hi_yi;
    }
    
    innovation = J_h_x * stateCovariance * J_h_x.transpose() /* +R */;
    kalmanGain = stateCovariance * J_h_x.transpose() * innovation.inverse();
    
    //compute aposteriori probability
    state = state + kalmanGain*y_k;
    stateCovariance = (Eigen::MatrixXd::Identity(stateCovariance.rows(),stateCovariance.cols())-kalmanGain*J_h_x)*stateCovariance;
    

}

void EigenEKFSLAM::addLandmark(MeasurementVector _stateFeature)
{
    double x0 = state[0], y0 = state[1], z0 = state[2], qw0 = state[3], qx0 = state[4], qy0 = state[5], qz0 = state[6];
    //landmakr absolute 3D point
    double xi = _stateFeature[0], yi = _stateFeature[1], zi = _stateFeature[2], qwi = _stateFeature[3], qxi = _stateFeature[4], qyi = _stateFeature[5], qzi = _stateFeature[6];
    Eigen::VectorXd newLandmark;
    newLandmark << x0 - xi*(2*pow(qy0,2) + 2*pow(qz0,2) - 1) - yi*(2*qw0*qz0 - 2*qx0*qy0) + zi*(2*qw0*qy0 + 2*qx0*qz0),
    y0 - yi*(2*pow(qx0,2) + 2*pow(qz0,2) - 1) + xi*(2*qw0*qz0 + 2*qx0*qy0) - zi*(2*qw0*qx0 - 2*qy0*qz0),
    z0 - zi*(2*pow(qx0,2) + 2*pow(qy0,2) - 1) - xi*(2*qw0*qy0 - 2*qx0*qz0) + yi*(2*qw0*qx0 + 2*qy0*qz0),
    qw0*qwi - qx0*qxi - qy0*qyi - qz0*qzi,
    qw0*qxi + qwi*qx0 + qy0*qzi - qyi*qz0,
    qw0*qyi + qwi*qy0 - qx0*qzi + qxi*qz0,
    qw0*qzi + qwi*qz0 + qx0*qyi - qxi*qy0;
    Eigen::VectorXd newState;
    newState << state, newLandmark;
    
    
    Eigen::MatrixXd J_yl_xv, J_yl_hl;
    J_yl_xv<<
    1, 0, 0, 2*qy0*zi - 2*qz0*yi,            2*qy0*yi + 2*qz0*zi, 2*qx0*yi - 4*qy0*xi + 2*qw0*zi, 2*qx0*zi - 2*qw0*yi - 4*qz0*xi,
    0, 1, 0, 2*qz0*xi - 2*qx0*zi, 2*qy0*xi - 4*qx0*yi - 2*qw0*zi,            2*qx0*xi + 2*qz0*zi, 2*qw0*xi - 4*qz0*yi + 2*qy0*zi,
    0, 0, 1, 2*qx0*yi - 2*qy0*xi, 2*qz0*xi + 2*qw0*yi - 4*qx0*zi, 2*qz0*yi - 2*qw0*xi - 4*qy0*zi,            2*qx0*xi + 2*qy0*yi,
    0, 0, 0,                 qwi,                           -qxi,                           -qyi,                           -qzi,
    0, 0, 0,                 qxi,                            qwi,                            qzi,                           -qyi,
    0, 0, 0,                 qyi,                           -qzi,                            qwi,                            qxi,
    0, 0, 0,                 qzi,                            qyi,                           -qxi,                            qwi;
    
    J_yl_hl<<
   -2*qy0*qy0 - 2*qz0*qz0 + 1,   2*qx0*qy0 - 2*qw0*qz0,   2*qw0*qy0 + 2*qx0*qz0,   0,    0,    0,    0,
    2*qw0*qz0 + 2*qx0*qy0, - 2*qx0*qx0 - 2*qz0*qz0 + 1,   2*qy0*qz0 - 2*qw0*qx0,   0,    0,    0,    0,
    2*qx0*qz0 - 2*qw0*qy0,   2*qw0*qx0 + 2*qy0*qz0, - 2*qx0*qx0 - 2*qy0*qy0 + 1,   0,    0,    0,    0,
                        0,                       0,                       0, qw0, -qx0, -qy0, -qz0,
                        0,                       0,                       0, qx0,  qw0, -qz0,  qy0,
                        0,                       0,                       0, qy0,  qz0,  qw0, -qx0,
                        0,                       0,                       0, qz0, -qy0,  qx0,  qw0;
    Eigen::MatrixXd newStateCov(newState.size(), newState.size());
    newStateCov.block(0,0, state.size(), state.size())= stateCovariance;
    for(int i = 0; i < state.size(); i++){
        newStateCov.block(state.size(), i*7, 7, 7) = J_yl_xv*stateCovariance.block(0,i*7,7,7);
        newStateCov.block(i*7, state.size(), 7, 7) = newStateCov.block(state.size(), i*7, 7, 7).transpose();
    }
    newStateCov.block(state.size(), state.size(), 7, 7) = J_yl_xv*stateCovariance.block(0,0,7,7)*J_yl_xv.transpose()/* + J_yl_hl*R*J_yl_hl.transpose()*/;
    
    state = newState;
    stateCovariance = newStateCov;
    
}
const Eigen::VectorXd EigenEKFSLAM::getState()
{
    return state;
}
