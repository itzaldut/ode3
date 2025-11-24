#include "RKn.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>

using namespace std;

struct Params {
    double g;   // gravity
    double m;   // mass
    double d;   // diameter
    double b;   // linear drag
    double c;   // quadratic drag
    double xend; // target distance
};

// RK4 functions for motion with air resistance
double f_x(double t, const vector<double> &y, void *params) { return y[1]; }  // dx/dt = vx
double f_z(double t, const vector<double> &y, void *params) { return y[3]; }  // dz/dt = vz

double f_vx(double t, const vector<double> &y, void *params) {
    Params *p = (Params*)params;
    double vx = y[1], vz = y[3];
    double v = sqrt(vx*vx + vz*vz);
    return -p->b * vx - p->c * v * vx; // drag along x
}

double f_vz(double t, const vector<double> &y, void *params) {
    Params *p = (Params*)params;
    double vx = y[1], vz = y[3];
    double v = sqrt(vx*vx + vz*vz);
    return -p->b * vz - p->c * v * vz - p->g; // drag + gravity
}

// Stop function: stop if z < 0
double f_stop(double t, const vector<double> &y, void *params){
    return (y[2] < 0) ? 1 : 0;
}

int main(int argc, char **argv){
    Params pars;
    pars.g = 9.81;
    pars.m = 0.145;
    pars.d = 0.075;
    pars.b = 1.6e-4 * pars.d;
    pars.c = 0.25 * pars.d * pars.d;
    pars.xend = 18.5;

    double z0 = 1.4;
    double theta0 = 1.0 * M_PI / 180.0; // radians

    // Optional command-line overrides
    int c;
    while((c = getopt(argc, argv, "x:z:t:")) != -1){
        switch(c){
            case 'x': pars.xend = atof(optarg); break;
            case 'z': z0 = atof(optarg); break;
            case 't': theta0 = atof(optarg) * M_PI / 180.0; break;
        }
    }

    // RK4 function pointers
    vector<pfunc_t> funcs(4);
    funcs[0] = f_x;
    funcs[1] = f_vx;
    funcs[2] = f_z;
    funcs[3] = f_vz;

    // Bisection to find vPitch
    double vmin = 0.0, vmax = 50.0, vPitch = 0.0;
    double tol = 1e-4;
    for(int iter=0; iter<50; iter++){
        double vtry = 0.5*(vmin + vmax);

        vector<double> y(4);
        y[0] = 0.0;          // x
        y[1] = vtry * cos(theta0); // vx
        y[2] = z0;           // z
        y[3] = vtry * sin(theta0); // vz

        double t0 = 0.0;
        auto traj = RK4SolveN(funcs, y, 2000, t0, 10.0, &pars, f_stop);

        double xfinal = y[0]; // after RK4 integration, y[0] is final x

        if(fabs(xfinal - pars.xend) < tol){
            vPitch = vtry;
            break;
        }
        if(xfinal < pars.xend) vmin = vtry;
        else vmax = vtry;
        vPitch = vtry; // last try if tolerance not exactly hit
    }

    printf("********************************\n");
    printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n", pars.xend, z0, theta0*180.0/M_PI);
    printf("v_pitch = %lf m/s\n", vPitch);
    printf("********************************\n");

    return 0;
}

