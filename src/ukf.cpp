#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  ///* time when the state is true, in us
  time_us_ = ;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* Weights of sigma points
  n_sigma_ = 2*n_aug_+1;
  //set weights
  for(int i=0;i<n_sigma_;i++){
    if(i==0){
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    }else{
      weights_(i) = .5 / (lambda_ + n_aug_);
    }
  }
  weights_ = ;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
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
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  //P_aug.bottomRightCorner((2,2)) = Q;
  P_aug(5,5) = pow(std_a_,2);
  P_aug(6,6) = pow(std_yawdd_,2);

  //calculate square root of P
  MatrixXd A_aug = P_.llt().matrixL();
  MatrixXd A_lambda = A_aug * sqrt(lambda_ + n_x_);

  Xsig_aug.col(0) = x_aug;
  for (int col=0; col<n_aug_; col++){
    Xsig_aug.col(col+1) = x_aug + A_lambda.col(col);
    Xsig_aug.col(col+1+n_x_) = x_aug - A_lambda.col(col);
  }

  VectorXd delta_x(n_x_);
  VectorXd noise(n_x_);

  //predict sigma points
  for(int col=0; col<n_sigma_; col++){
    float px = Xsig_aug(0,col);
    float py = Xsig_aug(1,col);
    float v = Xsig_aug(2,col);
    float yaw = Xsig_aug(3,col);
    float yaw_d = Xsig_aug(4,col);

    float delta_yaw = yaw_d*delta_t;

    if(yaw_d == 0){
      delta_x(0) = v * cos(yaw) * delta_t;
      delta_x(1) = v * sin(yaw) * delta_t;
      delta_x(2) = 0;
      delta_x(3) = delta_yaw;
      delta_x(4) = 0;
    }else{
      delta_x(0) = (v/yaw_d) * (sin(yaw + delta_yaw) - sin(yaw));
      delta_x(1) = (v/yaw_d) * (-cos(yaw + delta_yaw) + cos(yaw));
      delta_x(2) = 0;
      delta_x(3) = delta_yaw;
      delta_x(4) = 0;
    }

    float noise_a = Xsig_aug(5,col);
    float noise_yawdd = Xsig_aug(6,col);

    float half_pow_delta_t = (.5) * pow(delta_t,2);

    noise(0) = half_pow_delta_t * cos(yaw) * noise_a;
    noise(1) = half_pow_delta_t * sin(yaw) * noise_a;
    noise(2) = delta_t * noise_a;
    noise(3) = half_pow_delta_t * noise_yawdd;
    noise(4) = delta_t * noise_yawdd;

    Xsig_pred_.col(col) = Xsig_aug.col(col).head(n_x_) + delta_x + noise;



  x_.fill(0.0);
  P_.fill(0.0);

  for(int row=0; row<n_x_; row++){
    x_(row) = (Xsig_pred_.row(row)*weights_).sum();
  }

  for(int i=0;i<n_sigma_;i++){

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * (x_diff)*(x_diff).transpose();
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

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  for(int i=0; i<n_sigma_;i++){
    float px = Xsig_pred_(0,i);
    float py = Xsig_pred_(1,i);
    float v = Xsig_pred_.(2,i);
    float yawrate = Xsig_pred_(3,i);
    float yawrate_dot = Xsig_pred_(4,i);

    float rho = sqrt(pow(px,2)+pow(py,2));
    float phi = atan(py/px);
    float rho_dot = (px*cos(yawrate)*v + py*sin(yawrate)*v) / rho;

    Zsig.col(i) << rho, phi, rho_dot;
  }

  //transform sigma points into measurement space
  z_pred.fill(0.0);
  for(int row=0; row<n_z; row++){
    z_pred(row) = (Zsig.row(row)*weights_).sum();
  }
  //calculate mean predicted measurement

  S.fill(0.0);
  for(int i=0;i<n_sigma_;i++){

    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  pow(std_radr_,2),0,0,
        0,pow(std_radphi_,2),0,
        0,0,pow(std_radrd_,2);

  //calculate measurement covariance matrix S
  S += R;


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0;i<n_sigma_;i++){
    //Tc += weights.col(i)*(Xsig_pred.col(i) - x)*(Zsig.col(i) - z_pred).transpose();

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

  }

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_,n_z);
  K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();


}
