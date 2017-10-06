#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "ERROR: invalid data to calculate RMSE" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		rmse += (estimations[i] - ground_truth[i]).array().pow(2).matrix();
	}

	rmse = rmse.array() / estimations.size();
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	Hj << 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0;

	float px2py2 = powf(px, 2) + powf(py, 2);
	float sqrt_px2py2 = sqrtf(px2py2);
	float px2py2_32 = powf(px2py2, 3.0 / 2);

	if (fabs(px2py2)<0.000001) {
		std::cout << "ERROR: px and py = 0 while calculating Jacobian" << std::endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << px / sqrt_px2py2, py / sqrt_px2py2, 0, 0,
		-py / px2py2, px / px2py2, 0, 0,
		py*(vx*py - vy*px) / px2py2_32, px*(vx*py - vy*px) / px2py2_32, px / sqrt_px2py2, py / sqrt_px2py2;

	return Hj;
}
