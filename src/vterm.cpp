///
/// @file vterm_energy_terminal_fixed.cpp
/// @brief Projectile motion: RK4 without air + terminal velocity vs mass with air
/// @author ChatGPT
/// @date 24 Nov 2025
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include "TPaveText.h"

using namespace std;

struct Params {
    double g;
    double m;
    double air_k;
};

// --- RK4 functions ---
double f_ri(double x, const vector<double> &y, void *params=0){ (void)x; return y[1]; }
double f_rj(double x, const vector<double> &y, void *params=0){ (void)x; return y[3]; }

double f_vi(double x, const vector<double> &y, void *params=0){
    (void)x;
    Params *p = (Params*)params;
    if(p->air_k==0.0) return 0;
    return -p->air_k * sqrt(y[1]*y[1]+y[3]*y[3]) * y[1] / p->m;
}

double f_vj(double x, const vector<double> &y, void *params=0){
    (void)x;
    Params *p = (Params*)params;
    double air = (p->air_k==0.0) ? 0.0 : -p->air_k * sqrt(y[1]*y[1]+y[3]*y[3]) * y[3] / p->m;
    return air - p->g;
}

// Stop when projectile hits the ground
double f_stop_ground(double x, const vector<double> &y, void *params=0){ return (y[2] < 0) ? 1 : 0; }

// Terminal velocity simulation never stops early
double f_stop_terminal(double x, const vector<double> &y, void *params=0){ return 0; }

int main(int argc, char **argv){
    // --- Parameters ---
    Params pars;
    pars.g = 9.81;
    pars.air_k = 0.1;  
    void *p_par = (void*)&pars;

    double theta = 45;   // degrees
    double v0 = 100;     // m/s

    int c;
    while ((c = getopt(argc, argv, "v:t:m:k:")) != -1){
        switch(c){
            case 'v': v0 = atof(optarg); break;
            case 't': theta = atof(optarg); break;
            case 'm': pars.m = atof(optarg); break;
            case 'k': pars.air_k = atof(optarg); break;
            case '?': fprintf(stderr,"Unknown option `%c'.\n",optopt); break;
        }
    }

    TApplication theApp("App", &argc, argv);
    UInt_t dh = gClient->GetDisplayHeight()/2;  
    UInt_t dw = 1.1*dh;

    // --- RK4 function pointers ---
    vector<pfunc_t> v_fun(4);
    v_fun[0] = f_ri;
    v_fun[1] = f_vi;
    v_fun[2] = f_rj;
    v_fun[3] = f_vj;

    vector<double> y(4);
    y[0] = 0;
    y[1] = v0 * cos(theta*M_PI/180.0);
    y[2] = 0;
    y[3] = v0 * sin(theta*M_PI/180.0);

    // --- 1) No air resistance ---
    pars.air_k = 0.0;  
    pars.m = 10.0;     
    double x0 = 0.0;
    auto tgN_noair = RK4SolveN(v_fun, y, 200, x0, 20.0, p_par, f_stop_ground);

    // Energy calculation
    TGraph *gKinetic   = new TGraph();
    TGraph *gPotential = new TGraph();
    TGraph *gTotal     = new TGraph();
    int npoints = tgN_noair[0].GetN();
    for(int i=0;i<npoints;i++){
        double t,xpos,ypos,vx,vy;
        tgN_noair[0].GetPoint(i,t,xpos);
        tgN_noair[2].GetPoint(i,t,ypos);
        const Double_t* vxArr = tgN_noair[1].GetY();
        const Double_t* vyArr = tgN_noair[3].GetY();
        vx = vxArr[i];
        vy = vyArr[i];
        double kinetic = 0.5*pars.m*(vx*vx+vy*vy);
        double potential = pars.m*pars.g*ypos;
        gKinetic->SetPoint(i,t,kinetic);
        gPotential->SetPoint(i,t,potential);
        gTotal->SetPoint(i,t,kinetic+potential);
    }

    // --- Plot trajectory ---
    TCanvas *c2 = new TCanvas("c2","Trajectory no air",dw,dh);
    tgN_noair[2].Draw("AL*");
    c2->Draw();

    // --- Plot energies ---
    TCanvas *cEnergy = new TCanvas("cEnergy","Energy vs Time",dw,dh);
    gKinetic->SetLineColor(kBlue); gKinetic->SetMarkerColor(kBlue);
    gPotential->SetLineColor(kGreen+2); gPotential->SetMarkerColor(kGreen+2);
    gTotal->SetLineColor(kRed); gTotal->SetMarkerColor(kRed);
    gKinetic->Draw("AL*");
    gPotential->Draw("L* SAME");
    gTotal->Draw("L* SAME");
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(gKinetic,"Kinetic Energy","l");
    leg->AddEntry(gPotential,"Potential Energy","l");
    leg->AddEntry(gTotal,"Total Energy","l");
    leg->Draw();
    cEnergy->Draw();

    // --- 2) Terminal velocity vs mass ---
    TGraph *gVtermSim = new TGraph();
    TGraph *gVtermAnal = new TGraph();
    pars.air_k = 0.1;
    int idx = 0;

    for(double m=0.001; m<=10.0; m*=1.5){ // 1 g to 10 kg
        pars.m = m;
        y[0]=0; y[1]=0; y[2]=0; y[3]=0; 
        double x_air = 0.0;
        int nsteps = 5000;
        double xmax = 200.0;

        auto tgN_air = RK4SolveN(v_fun, y, nsteps, x_air, xmax, p_par, f_stop_terminal);

        int nPoints = tgN_air[3].GetN();
        if(nPoints < 10) continue; // skip invalid runs
        const Double_t* vY = tgN_air[3].GetY();
        int nLast = nPoints / 10;
        double sum = 0.0;
        for(int j=nPoints - nLast; j<nPoints; j++) sum += fabs(vY[j]);
        double vtermSimAvg = sum / nLast;

        gVtermSim->SetPoint(idx, m, vtermSimAvg);

        double vtermAnal = sqrt(m*pars.g/pars.air_k);
        gVtermAnal->SetPoint(idx, m, vtermAnal);

        idx++;
    }

    TCanvas *cVterm = new TCanvas("cVterm","Terminal velocity vs Mass",dw,dh);
    gVtermSim->SetTitle("Terminal velocity vs Mass;Mass [kg];Terminal velocity [m/s]");
    gVtermSim->SetLineColor(kMagenta); gVtermSim->SetMarkerColor(kMagenta);
    gVtermSim->Draw("AL*");

    gVtermAnal->SetLineColor(kBlue); gVtermAnal->SetLineWidth(2);
    gVtermAnal->Draw("L SAME");

    TLegend *leg2 = new TLegend(0.65,0.7,0.9,0.9);
    leg2->AddEntry(gVtermSim,"Simulated","p");
    leg2->AddEntry(gVtermAnal,"Analytical","l");
    leg2->Draw();
    cVterm->Draw();

    TCanvas *cDiscussion = new TCanvas("cDiscussion","Discussion",dw,dh);
    TPaveText *pt = new TPaveText(0.05,0.05,0.95,0.95,"NDC");
    pt->SetTextSize(0.03);
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->AddText("Discussion:");
    pt->AddText("Energy is conserved well, having a small number of samples is");
    pt->AddText("better because there is less distortion");
    pt->AddText("Compared with the analytical result to check for accuracy");
    pt->Draw();
    cDiscussion->Draw();

    TCanvas *cSave = new TCanvas("cSave","SavePDF");
    cSave->Print("../vterm.pdf["); // open multipage PDF
    cEnergy->Print("../vterm.pdf");
    cVterm->Print("../vterm.pdf");
    cDiscussion->Print("../vterm.pdf");
    cSave->Print("../vterm.pdf]"); // close multipage PDF
 
   // --- Save all graphs ---
    TFile *tf = new TFile("RKnDemo_fixed.root","recreate");
    for(unsigned i=0;i<v_fun.size();i++) tgN_noair[i].Write();
    gKinetic->Write("KineticEnergy_noAir");
    gPotential->Write("PotentialEnergy_noAir");
    gTotal->Write("TotalEnergy_noAir");
    gVtermSim->Write("TerminalVelocity_vs_Mass_Sim");
    gVtermAnal->Write("TerminalVelocity_vs_Mass_Analytical");
    tf->Close();

    cout << "Simulation complete.\n";
    theApp.SetIdleTimer(30,".q");
    theApp.Run();
}






















