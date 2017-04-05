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

  n_sigma_ = 2*n_aug_+1;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  ///* time when the state is true, in us
  time_us_ = 0.0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

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

  NIS_radar_ = 0.0;
  NIS_laser_ = 0.0;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  ///* Weights of sigma points
  weights_ = VectorXd(n_sigma_);
  for(int i=0;i<n_sigma_;i++){
    if(i==0){
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    }else{
      weights_(i) = .5 / (lambda_ + n_aug_);
    }
  }


}

UKF::~UKF() {}

void UKF::InitializeCovariance(){
  P_ <<    1,   0,    0,   0,   0,
           0,   1,    0,   0,   0,
           0,   0,    .1,  0,   0,
           0,   0,    0,   .1,  0,
           0,   0,    0,   0,   .1;

  return;
}

void UKF::InitialiseStateVector(MeasurementPackage meas_package){

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    float rho = meas_package.raw_measurements_[0];
    float phi = meas_package.raw_measurements_[1];
    float rho_dot = meas_package.raw_measurements_[2];

    float px = cos(phi) * rho;
    float py = sin(phi) * rho;
    float v= 0.0;
    float yaw = 0.0;
    float yaw_d = 0.0;
    x_ << px, py, v, yaw, yaw_d;

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    float px = meas_package.raw_measurements_[0];
    float py = meas_package.raw_measurements_[1];
    float v= 0.0;
    float yaw = 0.0;
    float yaw_d = 0.0;
    x_ << px, py, v, yaw, yaw_d;

  }
  return;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   *  Initialization
  ****************************************************************************/

  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;
    time_us_ = meas_package.timestamp_;


    //state covariance matrix P
    InitializeCovariance();
    InitialiseStateVector(meas_package);

    // done initializing, no need to predict or update
    if (x_[0] != 0 || x_[1] != 0) {
      is_initialized_ = true;
    }
    return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  //cout << dt << endl;
  time_us_ = meas_package.timestamp_;
  //cout << "Test Prediction" << endl;
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_){
        //cout << "Test Radar" << endl;
        UpdateRadar(meas_package);
      }
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    if (use_laser_){
      //cout << "Test Laser" << endl;
      UpdateLidar(meas_package);
    }
  }

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


  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

  GenerateAugmentedSigmaPoints(x_, P_, Xsig_aug);
  SigmaPointPrediction(delta_t, Xsig_aug, Xsig_pred_);
  PredictMeanAndCovariance(Xsig_pred_, x_, P_);

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
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = meas_package.raw_measurements_.rows();

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  PredictLidarMeasurement(Xsig_pred_, Zsig, z_pred, S);

  float px = meas_package.raw_measurements_[0];
  float py = meas_package.raw_measurements_[1];

  VectorXd z = VectorXd(n_z);
  z << px, py;

  VectorXd z_diff = UpdateState(z, z_pred, S, Xsig_pred_, Zsig, x_, P_);

  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

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
  int n_z = meas_package.raw_measurements_.rows();

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  PredictRadarMeasurement(Xsig_pred_, Zsig, z_pred, S);

  float rho = meas_package.raw_measurements_[0];
  float phi = meas_package.raw_measurements_[1];
  float rho_dot = meas_package.raw_measurements_[2];

  VectorXd z = VectorXd(n_z);
  z << rho, phi, rho_dot;

  VectorXd z_diff = UpdateState(z, z_pred, S, Xsig_pred_, Zsig, x_, P_);

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}

void UKF::GenerateAugmentedSigmaPoints(const VectorXd &x, const MatrixXd &P, MatrixXd &Xsig_aug) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(n_x_) = x;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P;
  //P_aug.bottomRightCorner((2,2)) = Q;
  P_aug(5,5) = pow(std_a_,2);
  P_aug(6,6) = pow(std_yawdd_,2);

  //calculate square root of P
  MatrixXd A_aug = P_aug.llt().matrixL();
  MatrixXd A_lambda = A_aug * sqrt(lambda_ + n_aug_);

  Xsig_aug.col(0) = x_aug;
  for (int col=0; col<n_aug_; col++){
    Xsig_aug.col(col+1) = x_aug + A_lambda.col(col);
    Xsig_aug.col(col+1+n_aug_) = x_aug - A_lambda.col(col);
  }
}

void UKF::SigmaPointPrediction(double delta_t, const MatrixXd &Xsig_aug, MatrixXd &Xsig_pred) {

  VectorXd delta_x(n_x_);
  VectorXd noise(n_x_);

  //predict sigma points
  for(int col=0; col<n_sigma_; col++) {
    float px = Xsig_aug(0, col);
    float py = Xsig_aug(1, col);
    float v = Xsig_aug(2, col);
    float yaw = Xsig_aug(3, col);
    float yaw_d = Xsig_aug(4, col);

    float delta_yaw = yaw_d * delta_t;

    if (fabs(yaw_d) < 0.001) {
      delta_x(0) = v * cos(yaw) * delta_t;
      delta_x(1) = v * sin(yaw) * delta_t;
      delta_x(2) = 0;
      delta_x(3) = delta_yaw;
      delta_x(4) = 0;
    } else {
      delta_x(0) = (v / yaw_d) * (sin(yaw + delta_yaw) - sin(yaw));
      delta_x(1) = (v / yaw_d) * (-cos(yaw + delta_yaw) + cos(yaw));
      delta_x(2) = 0;
      delta_x(3) = delta_yaw;
      delta_x(4) = 0;
    }

    float noise_a = Xsig_aug(5, col);
    float noise_yawdd = Xsig_aug(6, col);

    //cout << noise_yawdd << endl;
    //cout << noise_a << endl;

    float half_pow_delta_t = (.5) * pow(delta_t, 2);

    noise(0) = half_pow_delta_t * cos(yaw) * noise_a;
    noise(1) = half_pow_delta_t * sin(yaw) * noise_a;
    noise(2) = delta_t * noise_a;
    noise(3) = half_pow_delta_t * noise_yawdd;
    noise(4) = delta_t * noise_yawdd;

    Xsig_pred.col(col) = Xsig_aug.col(col).head(n_x_) + delta_x + noise;
  }
}

void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd &x, MatrixXd &P) {

  x.fill(0.0);
  for(int row=0; row<n_x_; row++){
    x(row) = (Xsig_pred.row(row)*weights_).sum();
  }

  P.fill(0.0);
  for(int i=0;i<n_sigma_;i++){
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
    x_diff(3) = normalizeRadiansPiToMinusPi(x_diff(3));
    //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P += weights_(i) * (x_diff)*(x_diff).transpose();
  }
}

void UKF::PredictLidarMeasurement(const MatrixXd &Xsig_pred, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S) {
    int n_z = z_pred.rows();

    //transform sigma points into measurement space
    for(int i=0; i<n_sigma_;i++){
        float px = Xsig_pred(0,i);
        float py = Xsig_pred(1,i);
        float v = Xsig_pred(2,i);
        float yawrate = Xsig_pred(3,i);
        float yawrate_dot = Xsig_pred(4,i);

        Zsig.col(i) << px, py;
    }

    //calculate mean predicted measurement
    z_pred.fill(0.0);
    for(int row=0; row<n_z; row++){
        z_pred(row) = (Zsig.row(row)*weights_).sum();
    }

    S.fill(0.0);
    for(int i=0;i<n_sigma_;i++){

        // state difference
        VectorXd z_diff = Zsig.col(i) - z_pred;

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    MatrixXd R = MatrixXd(n_z,n_z);
    R << pow(std_laspx_,2),0,
         0,pow(std_laspy_,2);

    //calculate measurement covariance matrix S
    S += R;
}

void UKF::PredictRadarMeasurement(const MatrixXd &Xsig_pred, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = z_pred.rows();

  for(int i=0; i<n_sigma_;i++){
    float px = Xsig_pred_(0,i);
    float py = Xsig_pred_(1,i);
    float v = Xsig_pred_(2,i);
    float yaw = Xsig_pred_(3,i);
    float yawrate = Xsig_pred_(4,i);

    float rho = sqrt(pow(px,2)+pow(py,2));
    float phi = atan(py/px);
    float rho_dot = 0.0;
    if (rho < 0.0001) {
      rho_dot = 0.0; //(px * cos(yaw) * v + py * sin(yaw) * v) / 0.0001;
    } else {
      rho_dot = (px * cos(yaw) * v + py * sin(yaw) * v) / rho;
    }
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
    z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R << pow(std_radr_,2),0,0,
       0,pow(std_radphi_,2),0,
       0,0,pow(std_radrd_,2);

  //calculate measurement covariance matrix S
  S += R;
}

auto UKF::UpdateState(const VectorXd &z, const VectorXd &z_pred, const MatrixXd &S,
                      const MatrixXd &Xsig_pred, const MatrixXd &Zsig, VectorXd &x, MatrixXd &P) -> VectorXd {

  int n_z = z_pred.rows();

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0;i<n_sigma_;i++){
    //Tc += weights.col(i)*(Xsig_pred.col(i) - x)*(Zsig.col(i) - z_pred).transpose();

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    x_diff(3) = normalizeRadiansPiToMinusPi(x_diff(3));
    //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

  }

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_,n_z);
  K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));

  //update state mean and covariance matrix
  x = x + K*z_diff;
  P = P - K*S*K.transpose();

  return z_diff;
}