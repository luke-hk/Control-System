#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <math.h>

using CppAD::AD;

size_t N = 10;
double dt = 0.1;

size_t x_start = 0;
size_t v_start = x_start + N;
size_t a_start = v_start + N;
size_t gasVol_start = a_start + N;
size_t thrust_start = gasVol_start + N - 1;

double CoD = .2;      //Coefficient of Drag
double csa = 1.0;     //Cross Sectional Area
double V_r = 1.0;     //Volume of Robot
double m = 1.0;       //Mass of Robot
double rho_f = 1.0;   //Density of Water
double setDepth = 1.0;

class FG_eval {
 public:
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(Eigen::VectorXd fg, const Eigen::VectorXd vars) {

    // Adding Cost section
    for (int t = 0; t < N; t++) {
      fg[0] += 10 * pow(vars[v_start + t], 2);
      fg[0] += 1500 * pow((vars[x_start + t] - setDepth), 2);
      fg[0] += pow(vars[a_start + t + 1] - vars[a_start + t], 2);  // Reference velocity cost
    }
  }
};