#!/usr/bin/env python
## @package RKnplotDemo
# Example for making plots from our differential equation solver.
#
# Requires that RKnDemo has been run first

import ROOT as r      # needed to load ROOT libraries
#import numpy as np    # only needed to use numpy arrays, see example below
import sys
from math import sqrt

# graphs are stored with generic names
# each of our 'n' dependent vars are plotted vs the independent variable (t)
tf=r.TFile("RKnDemo.root") # open file for read access
if tf.IsZombie():
    print("RKnDemo.root not found in current directory, run RKnDemo first")
    sys.exit(1)

tg_x_vs_t=tf.Get("xy0")   # dependent var[0] vs independent var
tg_vx_vs_t=tf.Get("xy1")  # etc...
tg_y_vs_t=tf.Get("xy2")
tg_vy_vs_t=tf.Get("xy3")

tc=r.TCanvas()
tc.Divide(2,2)

# it will often be intersting to instead plot one dependent var vs another
time=tg_x_vs_t.GetX()    # extract sample times
xval=tg_x_vs_t.GetY()    # extract array of positions along x axis
height=tg_y_vs_t.GetY()  # extract array of positions along y axis
vx=tg_vx_vs_t.GetY()     # extract x,y velocity
vy=tg_vy_vs_t.GetY()
nvals=tg_x_vs_t.GetN()

tc.cd(1)
tg_y_vs_x=r.TGraph(nvals,xval,height)   # make a new graph of y vs x!
tg_y_vs_x.SetTitle("Projectile y vs x;x [m];y [m]")
tg_y_vs_x.Draw("a*")

# plot energy of the system
m=1     # mass of projectile (must match definition in cpp!)
g=9.81  # gravity acceleration constant
tg_E_vs_t=r.TGraph()
tg_E_vs_t.SetTitle("Projectile Etot-E0 vs t;t [s];E [J]")
tg_v_vs_t=r.TGraph()
tg_v_vs_t.SetTitle("Projectile velocity vs t;t [s];v [m/s]")

e0=0    # initial energy of projectile
for i in range(nvals):
    v2=vx[i]*vx[i]+vy[i]*vy[i]
    etot=0.5*m*v2+m*g*height[i];
    if i==0: e0=etot
    tg_E_vs_t.SetPoint(i,time[i],etot-e0)
    tg_v_vs_t.SetPoint(i,time[i],sqrt(v2))

tc.cd(2)
tg_E_vs_t.Draw("a*")
tc.cd(3)
tg_v_vs_t.Draw("a*")

tc.Update()  # makes sure image on screen is up to date w/ all plot changes
tc.Print("Projectile.png") # or use Projectile.png, etc. to save your canvas


# note: if you would rather manipulate the data using something
# other than ROOT/TGraphs, you can convert the data into np arrays
# as follows:
# np_xval = np.fromiter(xval, dtype=np.float, count=nvals)
# np_height = np.fromiter(height, dtype=np.float, count=nvals)
# now you have nunmpy arrays to use for calculations/plots

print("Hit return to exit")
sys.stdin.readline()





