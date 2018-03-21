#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.75;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    is_initialized_ = false;
    
    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_aug_;
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    weights_ = VectorXd(n_x_, 2 * n_aug_ + 1);
    
    x_ = VectorXd (n_x_);
    x_ << 1, 1, 0, 0, 0.1;
   
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<   0.1, 0, 0, 0, 0,
            0, 0.1, 0, 0, 0,
            0, 0, 0.1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if (!is_initialized_) {
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            float rho = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];
            
            float px = rho * cos(phi);
            float py = rho * sin(phi);
            
            x_ << px, py, 0, 0, 0;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
            /**
             Initialize state.
             */
            float px = meas_package.raw_measurements_[0];
            float py = meas_package.raw_measurements_[1];
            
            x_ << px, py, 0, 0, 0;
        }
        
        double weight_0 = lambda_ / (lambda_ + n_aug_);
        weights_(0) = weight_0;
        
        for(int i=1; i<2*n_aug_+1; i++) {
            double weight_ = 0.5 / (lambda_ + n_aug_);
            weights_(i) = weight_;
        }
        
        time_us_ = meas_package.timestamp_;
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
    double delta_t_ = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    
    Prediction(delta_t_);
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        UpdateRadar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t_) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    
    //    Augmentation:
    VectorXd x_aug_ = VectorXd(n_aug_);
    MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
    MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    //create augmented mean state
    x_aug_.head(n_x_) = x_;
    x_aug_(5) = 0.0;
    x_aug_(6) = 0.0;
    
    //create augmented covariance matrix
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;
    P_aug_(5, 5) = std_a_ * std_a_;
    P_aug_(6, 6) = std_yawdd_ * std_yawdd_;
    
    //create square root matrix
    MatrixXd L_ = P_aug_.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug_.col(0) = x_aug_;
    for (int i=0; i<n_aug_; i++)
    {
        Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * L_.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L_.col(i);
    }
    
    //    Prediction:
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x_ = Xsig_aug_(0,i);
        double p_y_ = Xsig_aug_(1,i);
        double v_ = Xsig_aug_(2,i);
        double yaw_ = Xsig_aug_(3,i);
        double yawd_ = Xsig_aug_(4,i);
        double nu_a_ = Xsig_aug_(5,i);
        double nu_yawdd_ = Xsig_aug_(6,i);
        
        //predicted state values
        double px_p_, py_p_;
        
        //avoid division by zero
        if (fabs(yawd_) > 0.001) {
            px_p_ = p_x_ + v_ / yawd_ * (sin(yaw_ + yawd_ * delta_t_) - sin(yaw_));
            py_p_ = p_y_ + v_ / yawd_ * (cos(yaw_) - cos(yaw_ + yawd_ * delta_t_));
        }
        else {
            px_p_ = p_x_ + v_ * delta_t_ * cos(yaw_);
            py_p_ = p_y_ + v_ * delta_t_ * sin(yaw_);
        }
        
        double v_p_ = v_;
        double yaw_p_ = yaw_ + yawd_ * delta_t_;
        double yawd_p_ = yawd_;
        
        //add noise
        px_p_ = px_p_ + 0.5 * nu_a_ * delta_t_ * delta_t_ * cos(yaw_);
        py_p_ = py_p_ + 0.5 * nu_a_ * delta_t_ * delta_t_ * sin(yaw_);
        v_p_ = v_p_ + nu_a_ * delta_t_;
        
        yaw_p_ = yaw_p_ + 0.5 * nu_yawdd_ * delta_t_ * delta_t_;
        yawd_p_ = yawd_p_ + nu_yawdd_ * delta_t_;
        
        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p_;
        Xsig_pred_(1,i) = py_p_;
        Xsig_pred_(2,i) = v_p_;
        Xsig_pred_(3,i) = yaw_p_;
        Xsig_pred_(4,i) = yawd_p_;
    }
    
    //    Predicting Mean Covariance:
    //predict state mean
    x_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    
    //predict state covariance matrix
    P_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
        while(x_diff_(3) > M_PI) x_diff_(3) -= 2. * M_PI;
        while(x_diff_(3) < -M_PI) x_diff_(3) += 2. * M_PI;
        
        P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    int n_z_ = 2;
    
    VectorXd z_pred_ = VectorXd(n_z_);
    MatrixXd Zsig_ = MatrixXd(n_z_, 2*n_aug_+1);
    MatrixXd S_ = MatrixXd(n_z_, n_z_);
    MatrixXd R_ = MatrixXd(n_z_, n_z_);
    
    //create matrix for cross correlation Tc
    MatrixXd Tc_ = MatrixXd(n_x_, n_z_);
    
    //transform sigma points into measurement space
    for (int i=0; i<2*n_aug_+1; i++) {
        double p_x_ = Xsig_pred_(0, i);
        double p_y_ = Xsig_pred_(1, i);
        double v_ = Xsig_pred_(2, i);
        double yaw_ = Xsig_pred_(3, i);
        double v1_ = cos(yaw_) * v_;
        double v2_ = sin(yaw_) * v_;
        
        Zsig_(0, i) = sqrt(p_x_ * p_x_ + p_y_ * p_y_);
        Zsig_(1, i) = atan2(p_y_, p_x_);
        Zsig_(2, i) = (p_x_ * v1_ + p_y_ * v2_) / sqrt(p_x_ * p_x_ + p_y_ * p_y_);
    }
    
    //calculate mean predicted measurement
    z_pred_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
    }
    
    //calculate innovation covariance matrix S
    S_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
        
        while(z_diff_(1) > M_PI) z_diff_(1) -= 2. * M_PI;
        while(z_diff_(1) < -M_PI) z_diff_(1) += 2. * M_PI;
        
        S_ = S_ + weights_(i) * z_diff_ * z_diff_.transpose();
    }
    
    R_ <<  std_radr_ * std_radr_, 0, 0,
    0, std_radphi_ * std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;
    
    S_ = S_ + R_;
    
    //calculate cross correlation matrix
    Tc_.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
        //angle normalization
        while (z_diff_(1) > M_PI) z_diff_(1) -= 2. * M_PI;
        while (z_diff_(1) < -M_PI) z_diff_(1) += 2. * M_PI;
        
        // state difference
        VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        while (x_diff_(3) > M_PI) x_diff_(3) -= 2. * M_PI;
        while (x_diff_(3) < -M_PI) x_diff_(3) += 2. * M_PI;
        
        Tc_ = Tc_ + weights_(i) * x_diff_ * z_diff_.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K_ = Tc_ * S_.inverse();
    
    //residual
    VectorXd Z_ = meas_package.raw_measurements_;
    VectorXd z_diff_ = Z_ - z_pred_;
    
    //angle normalization
    while (z_diff_(1) > M_PI) z_diff_(1) -= 2. * M_PI;
    while (z_diff_(1) < -M_PI) z_diff_(1) += 2. * M_PI;
    
    //update state mean and covariance matrix
    x_ = x_ + K_ * z_diff_;
    P_ = P_ - K_ * S_ * K_.transpose();
    
    double nis_ = (Z_ - z_pred_).transpose() * S_.inverse() * (Z_ - z_pred_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    int n_z_ = 3;
    
    VectorXd z_pred_ = VectorXd(n_z_);
    MatrixXd Zsig_ = MatrixXd(n_z_, 2*n_aug_+1);
    MatrixXd S_ = MatrixXd(n_z_, n_z_);
    MatrixXd R_ = MatrixXd(n_z_, n_z_);
    
    //create matrix for cross correlation Tc
    MatrixXd Tc_ = MatrixXd(n_x_, n_z_);
    
    //transform sigma points into measurement space
    for (int i=0; i<2*n_aug_+1; i++) {
        double p_x_ = Xsig_pred_(0, i);
        double p_y_ = Xsig_pred_(1, i);
        double v_ = Xsig_pred_(2, i);
        double yaw_ = Xsig_pred_(3, i);
        double v1_ = cos(yaw_) * v_;
        double v2_ = sin(yaw_) * v_;
        
        if (p_x_ == 0 && p_y_ == 0) {
            Zsig_(0, i) = 0;
            Zsig_(1, i) = 0;
            Zsig_(2, i) = 0;
        }
        else {
            Zsig_(0, i) = sqrt(p_x_ * p_x_ + p_y_ * p_y_);
            Zsig_(1, i) = atan2(p_y_, p_x_);
            Zsig_(2, i) = (p_x_ * v1_ + p_y_ * v2_) / sqrt(p_x_ * p_x_ + p_y_ * p_y_);
        }
        
    }
    
    //calculate mean predicted measurement
    z_pred_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
    }
    
    //calculate innovation covariance matrix S
    S_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
        
        while(z_diff_(1) > M_PI) z_diff_(1) -= 2. * M_PI;
        while(z_diff_(1) < -M_PI) z_diff_(1) += 2. * M_PI;
        
        S_ = S_ + weights_(i) * z_diff_ * z_diff_.transpose();
    }
    
    R_ <<  std_radr_ * std_radr_, 0, 0,
    0, std_radphi_ * std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;
    
    S_ = S_ + R_;
    
    //calculate cross correlation matrix
    Tc_.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
        //angle normalization
        while (z_diff_(1) > M_PI) z_diff_(1) -= 2. * M_PI;
        while (z_diff_(1) < -M_PI) z_diff_(1) += 2. * M_PI;
        
        // state difference
        VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        while (x_diff_(3) > M_PI) x_diff_(3) -= 2. * M_PI;
        while (x_diff_(3) < -M_PI) x_diff_(3) += 2. * M_PI;
        
        Tc_ = Tc_ + weights_(i) * x_diff_ * z_diff_.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K_ = Tc_ * S_.inverse();
    
    //residual
    VectorXd Z_ = meas_package.raw_measurements_;
    VectorXd z_diff_ = Z_ - z_pred_;
    
    
    //angle normalization
    while (z_diff_(1) > M_PI) z_diff_(1) -= 2. * M_PI;
    while (z_diff_(1) < -M_PI) z_diff_(1) += 2. * M_PI;
    
    //update state mean and covariance matrix
    x_ = x_ + K_ * z_diff_;
    P_ = P_ - K_ * S_ * K_.transpose();
    
    double nis_ = (Z_ - z_pred_).transpose() * S_.inverse() * (Z_ - z_pred_);
}
