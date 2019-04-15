#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  if (estimations.size() == 0) {
    cout << "Estimations are empty." << endl;
    return rmse;
  }
  if (estimations.size() != ground_truth.size()) {
    cout << "Estimations and ground_truth should be of equal size." << endl;
    return rmse;
  }
  
  // accumulate squared residuals
  VectorXd sum(4);
  sum << 0, 0, 0, 0;
  
  for (int i=0; i < estimations.size(); ++i) {
    VectorXd resid(4);
    VectorXd sq_resid(4);
    resid = estimations[i] - ground_truth[i];
    
    sq_resid = resid.array()*resid.array();
    sum = sum + sq_resid;
  }
  
  // calculate the mean
  VectorXd mean(4);
  mean = sum / estimations.size();
  
  // calculate the squared root
  rmse = mean.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   *
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // check division by zero
  if ((px == 0) && (py == 0)) {
    cout << "CalculateJacobian() - Error - Division by Zero" << endl;
  }
  // compute the Jacobian matrix
  else {
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);
    Hj << (px/c2), (py/c2), 0, 0,
    -(py/c1), (px/c1), 0, 0,
    py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  }
  
  
  return Hj;
}
