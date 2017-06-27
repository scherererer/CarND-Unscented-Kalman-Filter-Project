#include "ukf.h"

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <cmath>
#include <iostream>
using std::vector;
using namespace std;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() :
	is_initialized_(false),
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_(true),
	// if this is false, radar measurements will be ignored (except during init)
	use_radar_(true),
	// initial state vector
	x_(5),
	// initial covariance matrix
	P_(5, 5),
	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_(30), ///< TODO This looks wildly off ....
	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_(30),
	// Laser measurement noise standard deviation position1 in m
	std_laspx_(0.15),
	// Laser measurement noise standard deviation position2 in m
	std_laspy_(0.15),
	// Radar measurement noise standard deviation radius in m
	std_radr_(0.3),
	// Radar measurement noise standard deviation angle in rad
	std_radphi_(0.03),
	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_(0.3),

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
	n_x_(5),
	n_aug_(7),
	weights_(2*n_aug_+1),
	lambda_(3 - n_aug_)
{
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (unsigned i = 1; i < weights_.rows (); ++i)
		weights_(i) = 1.0f / (2.0f * (lambda_ + n_aug_));
}

UKF::~UKF()
{
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage const &meas_package)
{
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/

	if (! is_initialized_)
	{
		Initialize (meas_package);
		return;
	}
}

void UKF::Initialize (MeasurementPackage const &meas_package)
{
	time_us_ = meas_package.timestamp_;

	switch (meas_package.sensor_type_)
	{
	case MeasurementPackage::LASER:
		// px
		// py
		// v
		// psi
		// psi-dot

		x_ << meas_package.raw_measurements_[0],
		      meas_package.raw_measurements_[1],
		      0,
		      0,
		      0;
		break;
	case MeasurementPackage::RADAR:
		float const p = meas_package.raw_measurements_[0];
		float const phi = meas_package.raw_measurements_[1];
		// rho-dot is ignored for initialization

		// px
		// py
		// v
		// psi
		// psi-dot

		x_ << p * cos (phi), // px
		      p * sin (phi), // py
		      0,
		      0,
		      0;
		break;
	}

	is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/

	/// \todo Generate sigma points
	MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

	Xsig.col (0) = x_;

	//calculate square root of P
	MatrixXd const A = P_.llt().matrixL();
	float const scale = sqrt (lambda_ + n_aug_);
	MatrixXd const sA = scale * A;

	for (unsigned i = 0; i < sA.cols (); ++i)
	{
		Xsig.col (i + 1) = x_ + sA.col (i);
		Xsig.col (i + 1 + sA.cols ()) = x_ - sA.col (i);
	}

	/// \todo Augment points
	//create augmented mean state
	VectorXd x_aug = VectorXd(7);

	x_aug.head (5) = x_;
	// mean of process noise is 0 so nothing to do

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	P_aug.topLeftCorner (5, 5) = P_;
	P_aug(5,5) = std_a_ * std_a_;
	P_aug(6,6) = std_yawdd_ * std_yawdd_;

	//create augmented sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	Xsig_aug.col (0) = x_aug;

	for (int i = 0; i < sA.cols (); ++i)
	{
		Xsig_aug.col (i + 1) = x_aug + sA.col (i);
		Xsig_aug.col (i + 1 + sA.cols ()) = x_aug - sA.col (i);
	}

	/// \todo Sigma point prediction

	/// \todo Predicted mean and covariance

	/// \todo Predict radar/lidar?
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage const &meas_package)
{
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
void UKF::UpdateRadar(MeasurementPackage const &meas_package)
{
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/
}
