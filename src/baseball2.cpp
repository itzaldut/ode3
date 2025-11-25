#include "RKn.hpp"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <unistd.h>
#include "TH1F.h"

using std::vector;
using std::cout;
using std::endl;

const double g = 32.174;           // ft/s^2
const double rho = 0.0023769;      // slugs/ft^3 (sea-level)
const double Cd = 0.35;            // drag coefficient (approx)
const double A = 0.00426;          // ft^2 (baseball cross-sectional area)
const double m = 0.32;             // slugs (~0.145 kg * 2.20462 / 32.174)
const double pi = 3.141592653589793;
const double S_MAG = 0.00041;   // chosen to get sensible break (ft units in force term)

//Parameters
struct Params {
    double wx, wy, wz;   // spin vector (rad/s)
};


double d_x(double t, const vector<double> &y, void *params) { return y[1]; }
double d_y(double t, const vector<double> &y, void *params) { return y[3]; }
double d_z(double t, const vector<double> &y, void *params) { return y[5]; }

double d_vx(double t, const vector<double> &y, void *params){
    Params *p = (Params*)params;
    double vx = y[1], vy = y[3], vz = y[5];
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v < 1e-8) return 0.0;
    double Fd = 0.5 * rho * v * v * Cd * A;
    double ux = vx / v;
    double magnus = S_MAG * (p->wy * vz - p->wz * vy);
    double ax = -(Fd/m) * ux + magnus;
    return ax;
}

double d_vy(double t, const vector<double> &y, void *params){
    Params *p = (Params*)params;
    double vx = y[1], vy = y[3], vz = y[5];
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v < 1e-8) return 0.0;
    double Fd = 0.5 * rho * v * v * Cd * A;
    double uy = vy / v;
    double magnus = S_MAG * (p->wz * vx - p->wx * vz);
    double ay = -(Fd/m) * uy + magnus;
    return ay;
}

double d_vz(double t, const vector<double> &y, void *params){
    Params *p = (Params*)params;
    double vx = y[1], vy = y[3], vz = y[5];
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v < 1e-8) return -g;
    double Fd = 0.5 * rho * v * v * Cd * A;
    double uz = vz / v;
    double magnus = S_MAG * (p->wx * vy - p->wy * vx);
    double az = -(Fd/m) * uz + magnus - g;
    return az;
}

//Initial Conditions
void SetupSlider(vector<double> &y0){
    double v0 = 125.0;             // ft/s (~85 mph)
    double theta = 1.0 * pi/180.0; 
    y0[0] = 0.0;          y0[1] = v0 * cos(theta);
    y0[2] = 0.0;          y0[3] = 0.0;
    y0[4] = 0.0;          y0[5] = v0 * sin(theta);
}

void SetupCurve(vector<double> &y0){
    double v0 = 125.0;             
    double theta = 1.0 * pi/180.0; 
    y0[0] = 0.0;          y0[1] = v0 * cos(theta);
    y0[2] = 0.0;          y0[3] = 0.0;
    y0[4] = 0.0;          y0[5] = v0 * sin(theta);
}

void SetupScrewball(vector<double> &y0){
    double v0 = 125.0;             // ft/s (~85 mph)
    double theta = 1.0 * pi/180.0;
    y0[0] = 0.0;          y0[1] = v0 * cos(theta);
    y0[2] = 0.0;          y0[3] = 0.0;
    y0[4] = 0.0;          y0[5] = v0 * sin(theta);
}

void SetupFastball(vector<double> &y0){
    double v0 = 140.0;             // ft/s (~95 mph)
    double theta = 1.0 * pi/180.0;
    y0[0] = 0.0;          y0[1] = v0 * cos(theta);
    y0[2] = 0.0;          y0[3] = 0.0;
    y0[4] = 0.0;          y0[5] = v0 * sin(theta);
}

// ---------- spin vector helper ----------
void SpinVector(double rpm, double phi_deg, Params &P){
    double w = rpm * 2.0 * pi / 60.0; // rad/s
    double phi = phi_deg * pi / 180.0;
    P.wx = 0.0;
    P.wy = w * sin(phi);
    P.wz = w * cos(phi);
}

int main(int argc, char **argv){
    bool showPlot = true;
    int c;
    while ((c = getopt(argc, argv, "n")) != -1) {
        if (c=='n') showPlot = false;
    }

    gStyle->SetOptStat(0);
    TApplication theApp("App",&argc,argv);

    const int Np = 4;
    const char *names[Np] = {"Slider","Curveball","Screwball","Fastball"};
    void (*setupFns[Np])(vector<double>&) = {SetupSlider, SetupCurve, SetupScrewball, SetupFastball};
    double rpms[Np] = {1800,1800,1800,1800};
    double phis[Np] = {0,45,135,225};
    int colors[Np] = {2,4,6,1};

    vector<pfunc_t> fnlist = {d_x, d_vx, d_y, d_vy, d_z, d_vz};

    TString pdfname = "../pitches.pdf";
    TCanvas c1("c1","Pitch trajectories",800,600);

    c1.Print(pdfname+"[");  // open PDF

    for(int ip=0; ip<Np; ++ip){
        vector<double> y(6,0.0);
        setupFns[ip](y);

        Params P;
        SpinVector(rpms[ip], phis[ip], P);

        vector<double> xs, ys, zs;
        xs.push_back(y[0]); ys.push_back(y[2]); zs.push_back(y[4]);

        double t=0.0, h=0.0005, tmax=5.0;
        while(y[0]<60.0 && t<tmax){
            y = RK4StepN(fnlist,y,t,h,&P);
            t += h;
            xs.push_back(y[0]); ys.push_back(y[2]); zs.push_back(y[4]);
        }

        double xlo=0, xhi=65, yabsmax=4, zhi=8;
        TH1F frame("frame","",100,xlo,xhi);
        frame.SetMinimum(-yabsmax); frame.SetMaximum(zhi);
        frame.SetXTitle("x (ft)"); frame.SetYTitle("y/z (ft)");
        frame.Draw();

        TGraph gxz(xs.size()), gxy(xs.size());
        for(size_t i=0;i<xs.size();++i){
            gxz.SetPoint(i,xs[i],zs[i]); // solid
            gxy.SetPoint(i,xs[i],ys[i]); // dashed
        }
        gxz.SetLineColor(colors[ip]); gxz.SetLineWidth(2); gxz.SetLineStyle(1);
        gxy.SetLineColor(colors[ip]); gxy.SetLineWidth(2); gxy.SetLineStyle(2);

        gxz.Draw("L SAME"); gxy.Draw("L SAME");

        TLegend leg(0.65,0.6,0.92,0.88);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(&gxz,names[ip],"l");
        leg.AddEntry(&gxy,Form("%s (y)",names[ip]),"l");
        leg.Draw();

        c1.Update();
        c1.Print(pdfname);  
    }

    c1.Print(pdfname+"]"); 

    if(showPlot){
        cout<<"Saved: "<<pdfname<<endl;
        theApp.Run();
    } else {
        cout<<"Saved: "<<pdfname<<" (no interactive display)"<<endl;
    }

    return 0;
}
