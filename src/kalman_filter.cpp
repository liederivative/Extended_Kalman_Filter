#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::_CorrectionStep(const VectorXd &y){
    
    MatrixXd Ht = H_.transpose();
//    std::cout<<"R_: "<<R_<<std::endl;
    MatrixXd PHt = P_ * Ht;
    MatrixXd S = H_ * PHt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
void KalmanFilter::Update(const VectorXd &z) {
    // predicted measurement
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    // Correction step
    this->_CorrectionStep(y);
    
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // calculate Jacobian Matrix
    
    // predicted measurement
    VectorXd z_pred = tools.Cartesian2Polar(x_);
    VectorXd y = z - z_pred;
    // check if between pi and -pi
    float phi = y(1);
    tools.Between_PI_PI(phi);
    y(1) = phi;
    // Correction Step
    this->_CorrectionStep(y);
}
