#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include <iostream>
#include <cmath>
#include <unistd.h>

using namespace std;

struct Params {
    double g;   // gravity
    double m;   // mass
    double d;   // diameter
    double b;   // air resistance
    double c;
};

// --- RK4 functions for 2D motion ---
double f_x(double x, const vector<double> &y, void *params=0){ return y[1]; }
double f_z(double x, const vector<double> &y, void *params=0){ return y[3]; }

double f_vx(double x, const vector<double> &y, void *params=0){
    Params *p = (Params*)params;
    double v = sqrt(y[1]*y[1]+y[3]*y[3]);
    return -p->b*v*y[1]/p->m;
}

double f_vz(double x, const vector<double> &y, void *params=0){
    Params *p = (Params*)params;
    double v = sqrt(y[1]*y[1]+y[3]*y[3]);
    return -p->b*v*y[3]/p->m - p->g;
}

// Stop when ball hits the ground
double f_stop(double x, const vector<double> &y, void *params=0){ return (y[2]<0) ? 1 : 0; }

double pitchHeightAtPlate(double v0, double theta0, double xend, Params &pars){
    double thetaRad = theta0*M_PI/180.0;
    vector<double> y(4);
    y[0] = 0;                // x position
    y[1] = v0*cos(thetaRad); // vx
    y[2] = pars.d;           // z position (release height)
    y[3] = v0*sin(thetaRad); // vz

    double x = 0;
    auto traj = RK4SolveN({f_x,f_vx,f_z,f_vz}, y, 500, x, 5.0, &pars, f_stop);
    return traj[2].GetY()[traj[2].GetN()-1]; // z at final x
}

int main(int argc, char **argv){
    Params pars;
    pars.g = 9.81;
    pars.m = 0.145;
    pars.d = 1.4; // example release height
    pars.b = 1.6e-4;
    pars.c = 0.25;

    double xend = 18.5;
    double theta0 = 1;  // degrees
    double zend_target = 0.9;

    // Bisection method to find vPitch
    double vLow = 10, vHigh = 50, vMid, zAtPlate;
    double tol = 1e-3;
    int maxIter = 50;

    for(int i=0; i<maxIter; i++){
        vMid = 0.5*(vLow+vHigh);
        zAtPlate = pitchHeightAtPlate(vMid, theta0, xend, pars);
        if(fabs(zAtPlate - zend_target)<tol) break;
        if(zAtPlate > zend_target) vHigh = vMid;
        else vLow = vMid;
    }

    printf("********************************\n");
    printf("(xend, z0, theta0) = (%lf, %lf, %lf)\n", xend, pars.d, theta0);
    printf("v_pitch = %lf m/s\n", vMid);
    printf("********************************\n");

    return 0;
}

