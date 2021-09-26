///
/// @file RKn.hpp 
/// @brief Runge-Kutta solvers. 
/// @author Bob Hirosky
/// @date 31 Dec 2019 
///
/// Generic solvers for linear first order ODEs using RK4.
///

// for convenience these functions return a graph of the solution for y(x)
#pragma once 
#include <TGraph.h>
#include <vector>
using std::vector;

// RK4 Solver for single ODE
TGraph RK4Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax);

// 
// RK4 solvers for set of N 1st order ODEs (1 independent parameter)
// The list of eqns is stored in a vector of fcn pointers
// initial conditions are given in a vector of doubles.


///
/// \var typedef double func_t(double, const vector<double>&)
/// \brief typedef double func_t(double, const vector<double>&)
///
/// B/c the 1st order ODEs may be coupled, each equation will have 
/// access to the entire set of dependent parameters, thus we use 
/// the interface: <br>
/// double func_t(double X, const vector<double> &Y),
/// \param[in] X: represents the independent paramater
/// \param[in] Y: is a vector of the dependent parameters
///
/// 
///
typedef double func_t(double, const vector<double>&, void *params);

///
/// \var typedef func_t* pfunc_t
/// \brief pointer to funct_t
///
typedef func_t* pfunc_t;

// The solvers will return a set of graphs for each dependent variable

TGraph RK4Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax);

vector<double> RK4StepN(vector<pfunc_t> &fnlist, vector<double> y,
			double x, double h, void *params=0);

vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
			 int nsteps, double x0, double xmax, void *params=0,
			 pfunc_t fstop=0);

void RK4SolveNx(vector<pfunc_t> &fnlist, vector<double> &y,
		int nsteps, double x0, double xmax, void *params=0,
		pfunc_t fstop=0);

vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
			 double h, double x0, pfunc_t fstop, void *params=0,
			 int nmax=1000);

vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
			  int nsteps, double x0, double xmax, void *params=0,
			  pfunc_t fstop=0, double errdef=1e-9, int maxrep=5);

vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
			  double h,  double x0, pfunc_t fstop, void *params=0,
			  double errdef=1e-9, int maxrep=5, int maxsteps=1000);

