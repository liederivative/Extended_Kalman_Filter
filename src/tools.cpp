#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
                    || estimations.size() == 0){
            std::cerr << "Invalid estimation or ground_truth data" << rmse << std::endl;
            exit(EXIT_FAILURE);
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){

            VectorXd residual = estimations[i] - ground_truth[i];

            //coefficient-wise multiplication
            residual = residual.array()*residual.array();
            rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
            std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
            c1 = 0.0001;
            //Hj.fill(0.0);
            
    }

    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
              -(py/c1), (px/c1), 0, 0,
              py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}
VectorXd Tools::Cartesian2Polar(const VectorXd& x_state){
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    float rho = hypot(px,py);
    
    float phi = atan2(py,px);
    
    if(fabs(rho) < 0.0001){
        std::cout << "Cartesian2Polar () - Error - Rho too small" << std::endl;
        rho = 0.0001;
    }
    
    float rho_dot = (px*vx+py*vy)/rho;
    
    VectorXd new_state(3);

    new_state << rho,phi,rho_dot;
//    std::cout << "new_state: " <<new_state << std::endl;
    return new_state;
}
VectorXd Tools::Polar2Cartesian(const VectorXd& polar_coords){
    //recover polar parameters
    float rho = polar_coords(0);
    float phi = polar_coords(1);
    float rho_dot = polar_coords(2);
    
    // calculate cartesian coords
    float px = rho*cos(phi);
    float py = rho*sin(phi);
    float vx = rho_dot * cos(phi);
    float vy = rho_dot * sin(phi);
    //return coords
    VectorXd new_state(4);
    new_state << px,py,vx,vy;
    return new_state;
}
void Tools::Between_PI_PI(float &angle){
  angle = atan2(sin(angle), cos(angle));
}