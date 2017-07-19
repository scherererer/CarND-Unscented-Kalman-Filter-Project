#include "ukf.h"

#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <cmath>
#include <iostream>
using namespace std;

namespace
{

float constexpr EPSILON = 0.001f;

inline float wrap (float a_)
{
	a_ = fmod (a_ + M_PI, 2.0 * M_PI);
	if (a_ < 0)
		a_ += M_PI * 2.0;
	return a_ - M_PI;
}

bool detectNan (string const &prefix_, MatrixXd const &A)
{
	bool hasnan = false;
	for (unsigned i = 0; i < A.rows (); ++i)
		for (unsigned j = 0; j < A.cols (); ++j)
			if (isnan (A(i,j)))
			{
				cerr << prefix_ << " NAN at index (" << i << "," << j << ")\n";
				hasnan = true;
			}
	return hasnan;
}

}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() :
	is_initialized_(false),
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_(true),
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
	std_a_(2),
	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_(1),
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
	if (! is_initialized_)
	{
		Initialize (meas_package);
		return;
	}

	if (! use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
		return;
	if (! use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
		return;

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
		cerr << "Skipping prediction due to 0s time difference\n";

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
	P_(2,2) = 15; // vel
	// Yaw can't possibly be off by more than M_PI.
	P_(3,3) = M_PI; // yaw
	// Can't imagine yaw rate being greater than M_PI rad/s
	P_(4,4) = M_PI; // yaw rate

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
	assert(delta_t > 0.0);

	//create augmented mean state
	VectorXd x_aug = VectorXd::Zero(7);

	assert(std::abs (x_(3)) <= M_PI);

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
	MatrixXd const A = P_aug.llt().matrixL();
	assert(! detectNan ("A", A));
	float const scale = sqrt(lambda_ + n_aug_);
	MatrixXd const sA = scale * A;

	MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

	Xsig_aug.col (0) = x_aug;

	for (int i = 0; i < sA.cols (); ++i)
	{
		Xsig_aug.col (i + 1) = x_aug + sA.col (i);
		Xsig_aug.col (i + 1 + n_aug_) = x_aug - sA.col (i);
	}

	assert(! detectNan ("Xsig_aug", Xsig_aug));

	//predict sigma points
	//avoid division by zero
	//write predicted sigma points into right column
	for (unsigned i = 0; i < Xsig_aug.cols (); ++i)
	{
		VectorXd const &x = Xsig_aug.col (i);

		float const px     = x(0);
		float const py     = x(1);
		float const v      = x(2);
		float const phi    = wrap (x(3));
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
		if (std::abs (phid) > EPSILON)
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
	}

	assert(! detectNan ("Xsig_pred_", Xsig_pred_));

	// Predicted mean
	x_.fill (0.0f);
	for (unsigned i = 0; i < weights_.rows (); ++i)
		x_ = x_ + (weights_(i) * Xsig_pred_.col(i));
	x_(3) = wrap(x_(3));

	/// Predict state covariance matrix
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
	int constexpr n_z = 2;

	VectorXd z = VectorXd (n_z);
	z << meas_package.raw_measurements_[0], // px
	     meas_package.raw_measurements_[1]; // py

	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	for (unsigned i = 0; i < Xsig_pred_.cols (); ++i)
	{
		float const px = Xsig_pred_.col(i)(0);
		float const py = Xsig_pred_.col(i)(1);

		// px
		Zsig.col(i)(0) = px;
		// py
		Zsig.col(i)(1) = py;
	}

	//calculate mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);

	for (unsigned i = 0; i < weights_.rows (); ++i)
		z_pred += (weights_(i) * Zsig.col(i));

	//calculate measurement covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);

	for (unsigned i = 0; i < weights_.rows (); ++i)
	{
		VectorXd zdiff = Zsig.col(i) - z_pred;

		S = S + (weights_(i) * zdiff * zdiff.transpose ());
	}

	MatrixXd R = MatrixXd(n_z, n_z);

	R << std_laspx_ * std_laspx_, 0,
	     0, std_laspy_ * std_laspy_;
	S = S + R;

	//calculate cross correlation matrix
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

	for (unsigned i = 0; i < weights_.rows (); ++i)
	{
		VectorXd zdiff = Zsig.col (i) - z_pred;
		VectorXd xdiff = Xsig_pred_.col (i) - x_;

		xdiff(3) = wrap(xdiff(3));

		Tc = Tc + (weights_ (i) * xdiff * zdiff.transpose ());
	}

	MatrixXd const Sinv = S.inverse ();

	//calculate Kalman gain K;
	MatrixXd const K = Tc * Sinv;
	//update state mean and covariance matrix
	VectorXd zdiff = z - z_pred;

	zdiff(1) = wrap(zdiff(1));
	x_ = x_ + K * zdiff;
	/// \note I'm doing this but it feels like I'm covering something else up
	x_(3) = wrap(x_(3));
	P_ = P_ - (K * S * K.transpose ());

	// Calculate NIS
	double const nis = zdiff.transpose () * Sinv * zdiff;
	cout << "LIDAR_NIS " << nis << "\n";
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage const &meas_package)
{
	int constexpr n_z = 3;

	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_[0], // p
		 wrap (meas_package.raw_measurements_[1]), // phi
		 meas_package.raw_measurements_[2]; // rho-dot


	for (unsigned i = 0; i < Xsig_pred_.cols (); ++i)
	{
		float const px = Xsig_pred_.col(i)(0);
		float const py = Xsig_pred_.col(i)(1);
		float const v = Xsig_pred_.col(i)(2);
		float const psi = Xsig_pred_.col(i)(3);
		float const phid = Xsig_pred_.col(i)(4);

		// rho
		Zsig.col(i)(0) = sqrt (px * px + py * py);
		// phi
		Zsig.col(i)(1) = atan2 (py, px);
		// rho-dot
		if (std::abs(Zsig.col(i)(0)) > EPSILON)
			Zsig.col(i)(2) = (px * cos(psi) * v + py * sin(psi) * v) / Zsig.col(i)(0);
		else
			// If rho is near 0 then assume it has no measurable radial velocity. I don't think
			// this should happen though because that would imply the target is on top of us.
			Zsig.col(i)(2) = 0.0f;
	}

	//calculate mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);

	for (unsigned i = 0; i < weights_.rows (); ++i)
		z_pred += (weights_(i) * Zsig.col(i));

	//calculate measurement covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);

	for (unsigned i = 0; i < weights_.rows (); ++i)
	{
		VectorXd zdiff = Zsig.col(i) - z_pred;

		zdiff(1) = wrap(zdiff(1));

		S = S + (weights_(i) * zdiff * zdiff.transpose ());
	}

	MatrixXd R = MatrixXd(n_z, n_z);

	R << std_radr_ * std_radr_, 0, 0,
	     0, std_radphi_ * std_radphi_, 0,
	     0, 0, std_radrd_ * std_radrd_;
	S = S + R;

	//calculate cross correlation matrix
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

	for (unsigned i = 0; i < weights_.rows (); ++i)
	{
		VectorXd zdiff = Zsig.col (i) - z_pred;

		zdiff(1) = wrap(zdiff(1));

		VectorXd xdiff = Xsig_pred_.col (i) - x_;

		xdiff(3) = wrap(xdiff(3));

		Tc = Tc + (weights_ (i) * xdiff * zdiff.transpose ());
	}

	MatrixXd const Sinv = S.inverse ();
	//calculate Kalman gain K;
	MatrixXd const K = Tc * Sinv;
	//update state mean and covariance matrix
	VectorXd zdiff = z - z_pred;

	zdiff(1) = wrap(zdiff(1));
	x_ = x_ + K * zdiff;
	/// \note I'm doing this but it feels like I'm covering something else up
	x_(3) = wrap(x_(3));
	P_ = P_ - (K * S * K.transpose ());

	// Calculate NIS
	double const nis = zdiff.transpose () * Sinv * zdiff;
	cout << "RADAR_NIS " << nis << "\n";
}
