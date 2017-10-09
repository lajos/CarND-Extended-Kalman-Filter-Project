#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "ERROR: invalid data to calculate RMSE" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {
		rmse += (estimations[i] - ground_truth[i]).array().pow(2).matrix();
	}

	rmse /= estimations.size();
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);

	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	Hj << 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0;

	double px2py2 = pow(px, 2) + pow(py, 2);
	double sqrt_px2py2 = sqrt(px2py2);
	double px2py2_32 = pow(px2py2, 3.0 / 2);

	if (abs(px2py2)<0.000001) {
		std::cout << "ERROR: px and py = 0 while calculating Jacobian" << std::endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << px / sqrt_px2py2, py / sqrt_px2py2, 0, 0,
		-py / px2py2, px / px2py2, 0, 0,
		py*(vx*py - vy*px) / px2py2_32, px*(vx*py - vy*px) / px2py2_32, px / sqrt_px2py2, py / sqrt_px2py2;

	return Hj;
}

VectorXd Tools::PolarToCartesian(const double& rho, const double& phi, const double& rhodot)
{
	double px = rho * cos(phi);
	double py = rho * sin(phi);
	double vx = rhodot * sin(phi);
	double vy = rhodot* cos(phi);

	VectorXd cartesian(4);
	cartesian << px, py, vx, vy;
	return cartesian;
}

VectorXd Tools::CartesianToPolar(const double& px, const double& py, const double& vx, const double& vy)
{
	double rho = sqrt(px * px + py * py);
	double phi = atan2(py, px);
	double rho_dot = (abs(rho) < 0.00001) ? 0 : (px * vx + py * vy) / rho;

	VectorXd polar(3);
	polar << rho, phi, rho_dot;
	return polar;
}

double Tools::ConstrainRadian(double x) {
	if (x > M_PI)
		return fmod(x - M_PI, 2 * M_PI) - M_PI;
	else if (x < -M_PI)
		return fmod(x + M_PI, 2 * M_PI) + M_PI;
	return x;
}
