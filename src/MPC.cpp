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
double g = 9.81;
//double setDepth = 1.0;

class FG_eval {
 public:
  Eigen::VectorXd coeffs;
  double setPoint;
  FG_eval(Eigen::VectorXd coeffs, double setPoint) { 
    this->coeffs = coeffs; 
    this->setPoint = setPoint;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    // Adding Cost section
    for (int t = 0; t < N; t++) {
      fg[0] += 5000 * CppAD::pow(vars[v_start + t], 2);
      fg[0] += 10000 * CppAD::pow((vars[x_start + t] - setPoint), 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);  // Reference velocity cost
      fg[0] += 500 * CppAD::pow((vars[thrust_start + t]), 2); // Limit thrust use
      fg[0] += 5000 * CppAD::pow((vars[gasVol_start + t + 1] - vars[gasVol_start + t]), 2); // Change in Gas Vol
    }

    // Initialization
    fg[1 + x_start] = vars[x_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + a_start] = vars[a_start];
    fg[1 + gasVol_start] = vars[gasVol_start];
    fg[1 + thrust_start] = vars[thrust_start];

    // Create predicted trajectory. Its length depends on N
    for (int t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> a1 = vars[a_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> gasVol0 = vars[gasVol_start + t - 1];
      AD<double> thrust0 = vars[thrust_start + t - 1];

      // Get next values
      double latency = 0.1;
      fg[1 + x_start + t] = x1 - (x0 + v0 * (dt + latency));
      fg[1 + v_start + t] = v1 - (v0 + a0 * (dt + latency));
      // hooplah = whole lotta;
      AD<double> PressDiff = x0/10 + 1;
      AD<double> B = 2 * rho_f * (V_r + gasVol0/PressDiff) * g;
      AD<double> Drag = CoD * csa * rho_f * pow(v0, 2);
      if (v0 >= 0) { //Moving down, so drag force is up
        fg[1 + a_start + t] = a1 - (thrust0 + (m * g) - B - Drag)/(2 * m); //Down is positive
      } else { //Moving up, so drag force is down
        fg[1 + a_start + t] = a1 - (thrust0 + (m * g) - B + Drag)/(2 * m); //Down is positive
      }
    }
  }
};

MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, double setPoint) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double v = state[1];
  double a = state[2];

  size_t n_vars = state.size() * N + 2 * (N-1);
  size_t n_constraints = N * state.size();

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  vars[x_start] = x;
  vars[v_start] = v;
  vars[a_start] = a;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (int i = 0; i < gasVol_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // Gas vol can't go less than 0 or greater than 50 mL
  for (int i = gasVol_start; i < thrust_start; i++) {
    vars_lowerbound[i] = 0;
    vars_upperbound[i] = 50;
  }

  // For thrust (4lbf = 17.5N)
  for (int i = thrust_start; i < n_vars; i++) {
    vars_lowerbound[i] = -17.5;
    vars_upperbound[i] = 17.5;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[a_start] = a;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[a_start] = a;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, setPoint);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  vector<double> result;
  result.push_back(solution.x[gasVol_start]);
  result.push_back(solution.x[thrust_start]);
  for (int i = 0; i < N-1; i++) {
    result.push_back(solution.x[x_start + i + 1]);
  }
  return result;
}