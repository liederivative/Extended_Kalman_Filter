#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

static const float PI = 3.14159265358979323844;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);
  /**
  * A helper method to convert cartesian coords to polar coords.
  */
  Eigen::VectorXd Cartesian2Polar(const Eigen::VectorXd& x_state);
  /**
  * A helper method to convert Polar coords to Cartesian coords.
  */
  Eigen::VectorXd Polar2Cartesian(const Eigen::VectorXd& polar_coords);
  /**
  * A helper method to keep angle between pi and -pi;
  */
  void Between_PI_PI( float &angle);
  
};

#endif /* TOOLS_H_ */
