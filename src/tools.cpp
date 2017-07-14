#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

#include <iostream>
using namespace std;


Tools::Tools() :
	rmse_ (4)
{
}

Tools::~Tools()
{
}

VectorXd Tools::CalculateRMSE(vector<VectorXd> const &estimations,
                              vector<VectorXd> const &ground_truth)
{
	/**
	TODO:
	* Calculate the RMSE here.
	*/

	rmse_.fill (0.0);

	if (estimations.empty () || ground_truth.empty ())
    {
        cerr << "Empty estimation or ground truth vector\n";
        return rmse_;
    }

    if (estimations.size () != ground_truth.size ())
    {
        cerr << "Estimation and ground truth vector must be same length\n";
        return rmse_;
    }

	// accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i)
	{
		VectorXd const diff = estimations[i] - ground_truth[i];
		VectorXd const squared = diff.array () * diff.array ();
		rmse_ = rmse_ + squared;
	}

	rmse_ = rmse_ / estimations.size ();

	return rmse_.array ().sqrt ();
}
