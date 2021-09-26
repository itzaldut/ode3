///
/// @file RKn.cpp 
/// @brief Runge-Kutta solvers. 
/// @author Bob Hirosky
///
/// Generic solvers for linear first order ODEs using RK4.
///

#include "RKn.hpp"
#include <assert.h>
#include <algorithm>
#include <cmath>

using std::max;
using std::min;

/// 
/// @fn TGraph RK4Solve(double (*f)(double x, double y), double y0,
///		int nsteps, double x0, double xmax);
/// @brief RK4 solver for a single ODE.
///
/// Solve an ODE of the form f = y'(x,y)
/// \param[in] f pointer to ODE function
/// the equation may depend on one independent variable (x) and a
/// dependent variable (y)
/// \param[in] y0 initial condition for dependent variable
/// \param[in] nsteps steps to use in the approximation
/// \param[in] x0 initial condition for independent variable
/// \param[in] xmax end of range for dependent variable in calculation
/// x0, xmax, nsteps  are used to set the step size
///
/// Returns a TGraph of y vs x.  x0 and y0 are not updated.
/// 
TGraph RK4Solve(double(*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax){

  double h=(xmax-x0)/nsteps;     // step size
  double x=x0;                   // independent variable
  
  double y=y0;          // dependent variable to plot vs x
  TGraph tg;
  tg.SetPoint(0,x0,y0);
	      
  double k1,k2,k3,k4;
  for (int i=0; i<nsteps-1; i++){
    k1 = h * f(x,y);
    k2 = h * f(x+h/2,y+k1/2);
    k3 = h * f(x+h/2,y+k2/2);
    k4 = h * f(x+h,y+k3);
    y = y + k1/6 + k2/3 + k3/3 + k4/6;
    x+=h;
    tg.SetPoint(i+1,x,y);
  }
  return tg;
}

/// @fn vector<double> RK4StepN(vector<pfunc_t> &fnlist, vector<double> y,
///			double x, double h, void* params)
/// @brief Perform one step using the RK4 algorithm
///
/// Call this function directly instead of using the higher level interfaces
/// for full control over the data in your intermediate steps.
/// @param[in] fnlist vector of function pointers to the ODEs describing the system
/// @param[in] y vector of initial conditions (this is updated by the step taken)
/// @param[in] x value of dependent vaaiable
/// @param[in] h step size to use
/// @param[in] params pointer to function parameters
///
/// Returns the new vector of y values
vector<double> RK4StepN(vector<pfunc_t> &fnlist, vector<double> y,
			double x, double h, void* params){
  int nFcn=fnlist.size();
  vector<double> k1(nFcn), k2(nFcn), k3(nFcn), k4(nFcn);
  vector<double> ytmp(nFcn);

  for (int i=0; i<nFcn; i++){
    k1[i] = h * fnlist[i](x, y, params);
    ytmp[i] = y[i] + k1[i]/2;
  }
  for (int i=0; i<nFcn; i++){
    k2[i] = h * fnlist[i](x+h/2, ytmp, params);
    ytmp[i] = y[i] + k2[i]/2;
  }
  for (int i=0; i<nFcn; i++){
    k3[i] = h * fnlist[i](x+h/2, ytmp, params);
    ytmp[i] = y[i] + k3[i];
  }
  for (int i=0; i<nFcn; i++){
    k4[i] = h * fnlist[i](x+h, ytmp, params);
  }
  // calculate next step
  for (int i=0; i<nFcn; i++) {
    y[i] = y[i] + k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6;
  }
  return y;
}

// make the graphs that will return information on the solution
vector<TGraph> SetupGraphs(int nFcn){
  vector<TGraph> tg(nFcn);
  for (int i=0; i<nFcn; i++){
    char buf[100];
    sprintf (buf,";independent variable (x);dependent variable y[%d]",i);
    tg[i].SetTitle(buf);
    sprintf(buf,"xy%d",i);
    tg[i].SetName(buf);
  }
  return tg;
}

///
/// \fn vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
///			 int nsteps, double x0, double xmax, void *params, pfunc_t fstop=0)
/// \brief RK4 solver for a system of ODEs.
///
/// Solves a series of coupled ODEs using the RK4 method using a fixed number of steps in range [x0,xmax]
/// \param[in] fnlist vector of function pointers to the ODEs describing the system
/// \param[in] y0 vector of initial conditions (need one per ODE)
/// \param[in] nsteps maximum number of steps in simulation
/// \param[in] x0 starting value of dependent vaaiable
/// \param[in] xmax maximum value of dependent variable
/// \param[in] params pointer to function parameters
/// \param[in] fstop optional function of dependent varibles to define a stopping
///        condition based on the results of the approximation
///
/// Returns a set of TGraphs of y_i vs x for each dependent variable
///
vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
			 int nsteps, double x0, double xmax, void *params,
			 pfunc_t fstop){

  assert(fnlist.size() == y0.size());   // need one initial condition per function
  int nFcn=fnlist.size();
  double h=(xmax-x0)/nsteps;     // step size
  double x=x0;                   // independent variable
  vector<double> y = y0;         // start with initial conditions on dependent vars
  vector<double> ytmp(nFcn);
  
  vector<TGraph> tg = SetupGraphs(nFcn);
  for (int i=0; i<nFcn; i++) tg[i].SetPoint(0,x0,y0[i]);
  
  for (int n=0; n<nsteps; n++){
    ytmp=RK4StepN(fnlist, y, x, h, params);
    // advance to next step and store results in graphs
    y=ytmp;
    x+=h;
    for (int i=0; i<nFcn; i++) tg[i].SetPoint(n+1,x,y[i]);
    if (fstop && fstop(x,ytmp,params)) break;
  }
  return tg;
}



///
/// \fn void RK4SolveNx(vector<pfunc_t> &fnlist, vector<double> &y,
///			 int nsteps, double x0, double xmax, void *params, pfunc_t fstop)
/// \brief RK4 solver for a system of ODEs.
///
/// Solves a series of coupled ODEs using the RK4 method using a fixed number of steps in range [x0,xmax]
/// \param[in] fnlist vector of function pointers to the ODEs describing the system
/// \param[in] y vector of initial conditions (updated in this function to the final location) 
/// \param[in] nsteps maximum number of steps in simulation
/// \param[in] x0 starting value of dependent variable
/// \param[in] xmax maximum value of dependent variable
/// \param[in] params pointer to function parameters
/// \param[in] fstop optional function of dependent varibles to define a stopping
///        condition based on the results of the approximation
///
///
void RK4SolveNx(vector<pfunc_t> &fnlist, vector<double> &y,
		int nsteps, double x0, double xmax, void *params, pfunc_t fstop){

  assert(fnlist.size() == y.size());   // need one initial condition per function
  int nFcn=fnlist.size();
  double h=(xmax-x0)/nsteps;     // step size
  double x=x0;                   // independent variable
  vector<double> ytmp(nFcn);
    
  for (int n=0; n<nsteps; n++){
    ytmp=RK4StepN(fnlist, y, x, h, params);
    y=ytmp;
    x+=h;
    //printf("%lg %lg %lg\n",x,y[0],y[1]);
    if (fstop && fstop(x,ytmp,params)) break;
  }
}


///
/// \fn vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
///			 double h, double x0, pfunc_t fstop, void *params, int nmax=1000);
/// \brief RK4 solver for a system of ODEs.
///
/// Solves a series of coupled ODEs using the RK4 method. This version is intended to run until the stopping condition is met.
/// \param[in] fnlist vector of function pointers to the ODEs describing the system
/// \param[in] y0 vector of initial conditions (need one per ODE)
/// \param[in] h step size to use
/// \param[in] x0 starting value of dependent vaaiable
/// \param[in] fstop function of dependent varibles to define a stopping
///        condition based on the results of the approximation
/// \param[in] params pointer to function parameters
/// \param[in] nmax maximum number of steps to use in case stopping isn't satisfied
///
/// Returns a set of TGraphs of y_i vs x for each dependent variable
///
vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
			 double h, double x0, pfunc_t fstop, void *params,
			 int nmax){
  double xmax=x0+h*nmax;  // serves as a backup stopping condition 
  return RK4SolveN(fnlist, y0, nmax, x0, xmax, params, fstop);
}

///
/// \fn vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
///			  int nsteps, double x0, double xmax, pfunc_t fstop=0,
///			  void *params, double errdef=1e-9, int maxrep=5);
/// \brief RK4 solver with adaptive step size for a system of ODEs.
///
/// \param[in] fnlist vector of function pointers to the ODEs describing the system
/// \param[in] y0 vector of initial conditions (need one per ODE)
/// \param[in] nsteps maximum number of steps in simulation
/// \param[in] x0 starting value of dependent vaaiable
/// \param[in] xmax maximum value of dependent variable
/// \param[in] params pointer to function parameters
/// \param[in] fstop function of dependent varibles to define a stopping
///        condition based on the results of the approximation
/// \param[in] errdef requested error per step. The error calculation
///         corresponds to a relative (absolute) error estimate
///         for y values large (small) compared to errdef
/// \param[in] maxrep max recalulations of h at each step to prevent infinite loops
///
/// Returns a set of TGraphs of y_i vs x for each dependent variable
///
/// This code loosely follows the examples in Fitzpatrick and
/// Numerical Recipes to adapt the step sizes
vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
			  int nsteps, double x0, double xmax, void *params, pfunc_t fstop,
			  double errdef, int maxrep){
  // **
  // The following could be passed as parameters, to allow more control
  // when calling the function but here they are hardcoded to (hopefully)
  // resonable values to simplify the function interface
  const double hmin=errdef*4;
  const double hmax=(xmax-x0)/4;
  // **

  assert(fnlist.size() == y0.size());   // need one initial condition per function
  int nFcn=fnlist.size();
  double h=(xmax-x0)/nsteps;     // step size
  double h_last=0, h_new=0;      // for tracking adapted step size
  double x=x0;                   // independent variable
  vector<double> y = y0;         // start w/ init. conditions on dependent vars
  vector<double> y1,y2;          // temporary storage for adapting steps
  
  vector<TGraph> tg = SetupGraphs(nFcn);
  for (int i=0; i<nFcn; i++) tg[i].SetPoint(0,x0,y0[i]);

  while (x<xmax){ // loop over steps in the x range
    int nreps=0;
    
    while (nreps<maxrep){ // take the steps, adapting step size as needed
      //full step
      y1=RK4StepN(fnlist, y, x, h, params);
      // two half steps
      y2=RK4StepN(fnlist, y, x, h/2, params);
      y2=RK4StepN(fnlist, y2, x+h/2, h/2, params);
      h_last=h;     // last step size used to evaluate the functions
      
      // Calculate truncation error, comparing the full step to two half steps
      double errAbs=0, errRel=0;
      for (int i=0; i<nFcn; i++){
	errAbs = max(errAbs,fabs(y1[i] - y2[i]));
	errRel = max(errRel,fabs(errAbs / y2[i]));
      }
      double err=max(errAbs,errRel);  // error estimate
      err=max(err,1e-15);             // protect against unlikely case of 0 error
      h_new = h * pow (fabs (errdef / err), 0.2);
      h_new=max(h_new,hmin);    // restrict adapted step size to above limits
      h_new=min(h_new,hmax);
      if (err<errdef) {
	h=h_new;
	break;                  // this step satisfies the requested errdef
      }
      nreps++;
      h=h_new;                  // if not try again
    }
    
    // advance to next step and store results in graphs
    y=y2;
    for (int i=0; i<nFcn; i++) tg[i].SetPoint(tg[i].GetN(),x+h_last,y[i]);
    if (fstop && fstop(x+h_last,y2,params)) break;  // fix me
    x+=h_new;
  }

  return tg;
}

///
/// \fn vector<TGraph>  RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
///			  double h,  double x0, pfunc_t fstop, void *params,
///			  double errdef=1e-9, int maxrep=5, int maxsteps=1000);
/// \brief RK4 solver with adaptive step size for a system of ODEs.
///
/// This version is intended to run until the stopping condition is met
/// \param[in] fnlist vector of function pointers to the ODEs describing the system
/// \param[in] y0 vector of initial conditions (need one per ODE)
/// \param[in] h initial guess for step size to use
/// \param[in] x0 starting value of dependent vaaiable
/// \param[in] fstop function of dependent varibles to define a stopping
///        condition based on the results of the approximation
/// \param[in] params pointer to function parameters
/// \param[in] errdef requested error per step. The error calculation
///         corresponds to a relative (absolute) error estimate
///         for y values large (small) compared to errdef
/// \param[in] maxrep max recalulations of h at each step to prevent infinite loops
/// \param[in] maxsteps maximum steps to use in approximation if stopping condition isn't satisfied 
///
/// Returns a set of TGraphs of y_i vs x for each dependent variable
///
vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
			  double h, double x0, pfunc_t fstop, void *params,
			  double errdef, int maxrep, int maxsteps){
  double xmax=x0+h*maxsteps;  // serves as a backup stopping condition
  return RK4SolveNA(fnlist, y0, maxsteps, x0, xmax, params, fstop, errdef, maxrep);
}
