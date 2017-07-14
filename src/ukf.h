#pragma once

#include "measurement_package.h"
#include "Eigen/Dense"

class UKF
{
public:
	/**
	* Constructor
	*/
	UKF();

	/**
	* Destructor
	*/
	virtual ~UKF();

	///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	Eigen::VectorXd const &x () const
		{ return x_; }

	/**
	* ProcessMeasurement
	* @param meas_package The latest measurement data of either radar or laser
	*/
	void ProcessMeasurement(MeasurementPackage const &meas_package);

private:
	/**
	 * Initialize the UKF
	 * @param meas_package Measurement package to use for initialization
	 */
	void Initialize (MeasurementPackage const &meas_package);

	/**
	* Prediction Predicts sigma points, the state, and the state covariance
	* matrix
	* @param delta_t Time between k and k+1 in s
	*/
	void Prediction(double delta_t);

	/**
	* Updates the state and the state covariance matrix using a laser measurement
	* @param meas_package The measurement at k+1
	*/
	void UpdateLidar(MeasurementPackage const &meas_package);

	/**
	* Updates the state and the state covariance matrix using a radar measurement
	* @param meas_package The measurement at k+1
	*/
	void UpdateRadar(MeasurementPackage const &meas_package);

	/// \brief Dump state of x_ for debugging
	void dumpState (std::string const &prefix_);

	///* initially set to false, set to true in first call of ProcessMeasurement
	bool is_initialized_;

	///* if this is false, laser measurements will be ignored (except for init)
	bool use_laser_;

	///* if this is false, radar measurements will be ignored (except for init)
	bool use_radar_;

	///* State dimension
	int const n_x_;

	///* Augmented state dimension
	int const n_aug_;

	///* Weights of sigma points
	Eigen::VectorXd weights_;

	///* Sigma point spreading parameter
	double const lambda_;

	///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	Eigen::VectorXd x_;

	///* state covariance matrix
	Eigen::MatrixXd P_;

	///* predicted sigma points matrix
	Eigen::MatrixXd Xsig_pred_;

	///* time when the state is true, in us
	long long time_us_;

	///* Process noise standard deviation longitudinal acceleration in m/s^2
	double const std_a_;

	///* Process noise standard deviation yaw acceleration in rad/s^2
	double const std_yawdd_;

	///* Laser measurement noise standard deviation position1 in m
	double const std_laspx_;

	///* Laser measurement noise standard deviation position2 in m
	double const std_laspy_;

	///* Radar measurement noise standard deviation radius in m
	double const std_radr_;

	///* Radar measurement noise standard deviation angle in rad
	double const std_radphi_;

	///* Radar measurement noise standard deviation radius change in m/s
	double const std_radrd_ ;
};
