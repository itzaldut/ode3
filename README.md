# ode3

physx630 odelib
---

To build the ODE library and example programs, simply type `make` in this top level ode3 directory.

Description of example programs:<br>

**RKnTest**: Solves a single 1st order ODE using the single equation RK4 solver and the ODE array solver

**RKnStep**: A basic example of the ODE array solver is applied to projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>. At each step in the elapsed time and x,y positions are printed.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m

Results for when I changed the values for baseball 1

**RKnDemo**: Solves for projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>
This program includes graphical output.  Detailed output is saved in TGraph objects in RKnDemo.root.  The file **RKnPlotDemo.py** shows how to access date in the TGraphs and can be used to generate additional plots.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m

**baseball1**:  Starter template for first baseball problem
(xend,z0,theta0) = (18.500000,1.400000,1.000000)
v_pitch = 31.369452 m/s
********************************
(my_root_env) [mgf3tf src]$ ./baseball1 -x 20 -t 45
********************************
(xend,z0,theta0) = (20.000000,1.400000,45.000000)
v_pitch = 13.689232 m/s
********************************
(my_root_env) [mgf3tf src]$ ./baseball1 -z 5
********************************
(xend,z0,theta0) = (18.500000,5.000000,1.000000)
v_pitch = 17.858696 m/s

**baseball2**:  Starter template for second baseball problem

**baseball_drag.ipynb**: this notebook describes the drag force equations used in the text.

gsl starter code
---

The starter code here (projGSL.cpp) demonstrates very basic usage of the gsl for solving a problem of coupled differential equations. An 8th order R-K solver with fixed step size is used. You are encouraged to try other solvers as you explore the problem. See here for the gsl docs: https://www.gnu.org/software/gsl/doc/html/ode-initval.html

This example solves the 2D projectile motion problem with a simple model for air resistance. After each step, data are stored in ROOT TGraphs, which are then displayed at the conclusion of the calculation.

The gsl provides a number of ODE solvers and a variety of interfaces.  Some of the solvers (not R-K methods) use the Jacobian matrix, which gives the devivative of the function wrt the dependent parameters.  See the gsl examples for details.

Python starter code
---

Two examples are given for using ODE solvers from the scipy.integrate sub-package in Python. In these examples graphs are made using matplotlib.

    Solution (projScPY2.py[ipynb]) using a more modern interface scipy.integrate.solve_ivp. See also: https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html and https://www.programcreek.com/python/example/119375/scipy.integrate.solve_ivp

    Solution (projScPY.py[ipynb]) using an older interface scipy.integrate.odeint¶ (see comments here: https://docs.scipy.org/doc/scipy/reference/integrate.html).  I do not recommend using this interface any longer.

The notebook versions contain additional comments on using the integrators.
