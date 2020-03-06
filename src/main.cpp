#include <math.h>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"

using Eigen::VectorXd;

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
}

int main() {

    // MPC is initialized here!
    MPC mpc;

    double x = 0.0;
    double v = 0.0;
    double a = 0.0;

    double setPoint = 1.5;

    int waypoints_number = 10;

    VectorXd time(waypoints_number);
    VectorXd waypoints_x(waypoints_number);

    auto coeffs = polyfit(time, waypoints_x, 3);

    Eigen::VectorXd state(3);
    // Current state
    state << x, v, a;

    auto vars = mpc.Solve(state, coeffs, setPoint);
    double gasVol = vars[0];
    double thrust = vars[1];

    // MPC Predicted trajectory
    vector<double> mpc_x_vals;

    for (int i = 2; i < vars.size(); i++) {
        mpc_x_vals.push_back(vars[i]);
    }

}