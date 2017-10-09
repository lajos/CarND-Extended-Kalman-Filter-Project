#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // laser measurement covariance matrix
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // radar measurement covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // laser state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;

  // laser measurement matrix
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  // radar Jacobian matrix
  Hj_ = MatrixXd(3, 4);

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  //initialize Q_
  ekf_.Q_ = MatrixXd(4, 4);

  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 1, 1, 1, 1;

  // initialize identity matrix
  long x_size = ekf_.x_.size();
  ekf_.I_ = MatrixXd::Identity(x_size, x_size);


  //set the acceleration noise components
  noise_ax = 5;
  noise_ay = 5;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // convert radar from polar to cartesian coordinates and initialize state.
		ekf_.x_ << tools.PolarToCartesian(measurement_pack.raw_measurements_[0],   // rho
			                               measurement_pack.raw_measurements_[1],   // phi
			                               measurement_pack.raw_measurements_[2]);  // rhodot
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		ekf_.x_ << measurement_pack.raw_measurements_[0],   // px
			       measurement_pack.raw_measurements_[1],   // py
			       0,                                       // vx (unknown)
			       0;                                       // vy (unknown)
	}

	// initialize timestamp
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // calculate elapsed time in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // update F_ matrix with elapsed time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // set process covariance matrix
  double dt44 = pow(dt, 4) / 4;
  double dt32 = pow(dt, 3) / 2;
  double dt2 = pow(dt, 2);
  ekf_.Q_ << dt44*noise_ax, 0, dt32*noise_ax, 0,
	  0, dt44*noise_ay, 0, dt32*noise_ay,
	  dt32*noise_ax, 0, dt2*noise_ax, 0,
	  0, dt32*noise_ay, 0, dt2*noise_ay;

  // Kalman predict
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // radar update
	  ekf_.H_ = Hj_;
	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
	  // laser update
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
