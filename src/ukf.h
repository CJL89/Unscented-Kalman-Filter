#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* Process noise covariance matrix for augmented state covariance matrix expand
  MatrixXd Q_;

  ///* Lidar Measurement noise covariance matrix
  MatrixXd R_lidar_;

  ///* Radar Measurement noise covariance matrix
  MatrixXd R_radar_;

  ///* predicted process sigma points matrix
  MatrixXd Xsig_pred_;
  ///* augmented state matrix
  MatrixXd Xsig_aug_;

  ///* time when the state is true, in us
  long long previous_time;
  double delta_t;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Radar sigma point in measurement space
  MatrixXd Zsig_radar_;

  ///* Lidar sigma point in measurement space
  MatrixXd Zsig_lidar_;

  ///* Measured mean states for radar
  VectorXd z_pred_radar_;

  ///* Measured mean states for lidar
  VectorXd z_pred_lidar_;

  ///* Radar innovation covariance matrix
  MatrixXd S_radar_;

  ///* Lidar innovation covariance matrix
  MatrixXd S_lidar_;

  ///* Radar measurement dimension, radar can measure rho, phi, and rho_dot
  int n_z_radar_;

  ///* Lidar measurement dimension, lidar measures px, py
  int n_z_lidar_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance matrix
   */
  void Prediction();

  /**
   * Calculate predicted mean and covariance of states
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Perform the update of state and covariance of Radar measurement
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
