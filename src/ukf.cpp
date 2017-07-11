#include "ukf.h"

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <cmath>
#include <iostream>
using std::vector;
using namespace std;

namespace
{

float const EPSILON = 0.001f;

inline float wrap (float a_)
{
	a_ = fmod (a_ + M_PI, 2.0 * M_PI);
	if (a_ < 0)
		a_ += M_PI * 2.0;
	return a_ - M_PI;
}

}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() :
	is_initialized_(false),
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_(false),
	// if this is false, radar measurements will be ignored (except during init)
	use_radar_(true),
	n_x_(5),
	n_aug_(7),
	weights_(2 * n_aug_ + 1),
	lambda_(3 - n_aug_),
	// initial state vector
	x_(5),
	// initial covariance matrix
	P_(Eigen::MatrixXd::Identity (5, 5)),
	// predicted sigma points matrix
	Xsig_pred_(MatrixXd::Zero(n_x_, 2 * n_aug_ + 1)),
	time_us_(0),
	// Process noise standard deviation longitudinal acceleration in m/s^2
	//std_a_(30), ///< TODO This looks wildly off ....
	std_a_(30), ///< TODO This is a random guess...
	// Process noise standard deviation yaw acceleration in rad/s^2
	//std_yawdd_(30),
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
	std_radrd_(0.3)

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
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

	if (! use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
		return;
	if (! use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
		return;

	cout << meas_package.timestamp_ / 1000000.0 << " Processing Measurement..." << endl;

	double const dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;

	/// \todo Handle going backwards in time more gracefully
	if (dt < 0.0)
	{
		cerr << "Warning: Time moved backwards by " << dt << " seconds. Resetting filter.\n";
		Initialize (meas_package);
		return;
	}

	if (dt > 0.0)
		Prediction (dt);
	else
		cout << "Skipping prediction due to 0s time difference" << endl;

	switch (meas_package.sensor_type_)
	{
	case MeasurementPackage::LASER:
		UpdateLidar (meas_package);
		break;
	case MeasurementPackage::RADAR:
		UpdateRadar (meas_package);
		break;
	}
}

void UKF::Initialize (MeasurementPackage const &meas_package)
{
	time_us_ = meas_package.timestamp_;
	P_.fill(0.0f);
	P_(0,0) = 1; // x
	P_(1,1) = 1; // y
	P_(2,2) = 1000; // vel
	P_(3,3) = 1000; // yaw
	P_(4,4) = 1000; // yaw rate

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
	assert (delta_t > 0.0);

	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/

	/// \todo Generate augmented sigma points
	//create augmented mean state
	VectorXd x_aug = VectorXd::Zero(7);

	x_aug.head(5) = x_;
	x_aug(5) = 0.0f;
	x_aug(6) = 0.0f;
	// mean of process noise is 0 so nothing to do

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd::Zero(7, 7);

	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5,5) = std_a_ * std_a_;
	P_aug(6,6) = std_yawdd_ * std_yawdd_;

	//create augmented sigma points
	//calculate square root of P
	/// \todo NaN risk here?
	MatrixXd const A = P_aug.llt().matrixL();
	float const scale = sqrt(lambda_ + n_aug_);
	MatrixXd const sA = scale * A;

	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	Xsig_aug.col (0) = x_aug;

	for (int i = 0; i < sA.cols (); ++i)
	{
		Xsig_aug.col (i + 1) = x_aug + sA.col (i);
		Xsig_aug.col (i + 1 + n_aug_) = x_aug - sA.col (i);
	}

	//predict sigma points
	//avoid division by zero
	//write predicted sigma points into right column
	for (unsigned i = 0; i < Xsig_aug.cols (); ++i)
	{
		VectorXd const &x = Xsig_aug.col (i);

		float const px     = x(0);
		float const py     = x(1);
		float const v      = x(2);
		float const phi    = x(3);
		float const phid   = x(4);
		float const va     = x(5);
		float const vphidd = x(6);

		VectorXd noise = VectorXd (5);
		noise << 0.5f * delta_t * delta_t * cos (phi) * va,
				 0.5f * delta_t * delta_t * sin (phi) * va,
				 delta_t * va,
				 0.5f * delta_t * delta_t * vphidd,
				 delta_t * vphidd;

		// if phi dot is not zero
		if (fabs (phid) > EPSILON)
		{
			VectorXd t2 = VectorXd (5);
			t2 << (v / phid) * (sin (phi + phid * delta_t) - sin (phi)),
			      (v / phid) * (-cos (phi + phid * delta_t) + cos (phi)),
			      0,
			      phid * delta_t,
			      0;

			Xsig_pred_.col (i) = x.head (5) + t2 + noise;
		}
		else
		{
			VectorXd t2 = VectorXd (5);
			t2 << v * cos (phi) * delta_t,
			      v * sin (phi) * delta_t,
			      0.0f,
			      phid * delta_t,
			      0.0f;

			Xsig_pred_.col (i) = x.head (5) + t2 + noise;
		}

		/// \todo Not sure if this is necessary
		//Xsig_pred_.col(i)(3) = wrap(Xsig_pred_.col(i)(3));
	}

	/// \todo Predicted mean and covariance
	x_.fill (0.0f);
	for (unsigned i = 0; i < weights_.rows (); ++i)
		x_ = x_ + (weights_(i) * Xsig_pred_.col(i));
	/// \todo not sure if this wrapping is necessary
	//x_(3) = wrap(x_(3));
	//predict state covariance matrix
	P_.fill (0.0f);
	for (unsigned i = 0; i < weights_.rows (); ++i)
	{
		VectorXd xdiff = Xsig_pred_.col (i) - x_;

		xdiff(3) = wrap(xdiff(3));

		P_ = P_ + (weights_(i) * xdiff * xdiff.transpose ());
	}
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
	cout << "UPDATE LIDAR" << endl;
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

	int constexpr n_z = 3;

	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	VectorXd z_pred = VectorXd::Zero(n_z);
	MatrixXd S = MatrixXd::Zero(n_z,n_z);
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_[0], // p
		 wrap (meas_package.raw_measurements_[1]), // phi
		 meas_package.raw_measurements_[2]; // rho-dot


	for (unsigned i = 0; i < Xsig_pred_.cols (); ++i)
	{
		float const px = Xsig_pred_.col(i)(0);
		float const py = Xsig_pred_.col(i)(1);
		float const v = Xsig_pred_.col(i)(2);
		float const phi = Xsig_pred_.col(i)(3);
		float const phid = Xsig_pred_.col(i)(4);

		// rho
		Zsig.col(i)(0) = sqrt (px * px + py * py);
		// phi
		Zsig.col(i)(1) = wrap (atan2 (py, px));
		// rho-dot
		/// \todo divide by zero risk?
		Zsig.col(i)(2) = (px * cos (phi) * v + py * sin (phi) * v) / Zsig.col(i)(0);
	}

	//calculate mean predicted measurement
	for (unsigned i = 0; i < weights_.rows (); ++i)
		z_pred = z_pred + (weights_(i) * Zsig.col(i));
	z_pred(1) = wrap (z_pred (1));

	//calculate measurement covariance matrix S
	MatrixXd R = MatrixXd::Zero(3, 3);

	R << std_radr_ * std_radr_, 0, 0,
	     0, std_radphi_ * std_radphi_, 0,
	     0, 0, std_radrd_ * std_radrd_;
	S = S + R;

	for (unsigned i = 0; i < weights_.rows (); ++i)
	{
		VectorXd zdiff = Zsig.col(i) - z_pred;

		zdiff(1) = wrap(zdiff(1));

		S = S + (weights_(i) * zdiff * zdiff.transpose ());
	}

	//calculate cross correlation matrix
	for (unsigned i = 0; i < weights_.rows (); ++i)
	{
		VectorXd zdiff = Zsig.col (i) - z_pred;

		zdiff(1) = wrap(zdiff(1));

		VectorXd xdiff = Xsig_pred_.col (i) - x_;

		xdiff(3) = wrap(xdiff(3));

		Tc = Tc + (weights_ (i) * xdiff * zdiff.transpose ());
	}
	//calculate Kalman gain K;
	MatrixXd const K = Tc * S.inverse ();
	//update state mean and covariance matrix
	VectorXd zdiff = z - z_pred;

	zdiff(1) = wrap(zdiff(1));
	x_ = x_ + K * zdiff;
	P_ = P_ - (K * S * K.transpose ());
}
