"""
Title: FD Solver for ML-2DOF Oscillator with friction
Author: Sebastian McDonald
Date: Feb 29, 2024    
    
Description: FD solver for pure friction nonlinearity of 2DOF oscillator
Reference: Mathias Legrand, Christophe Pierre. A compact, equality-based weighted
residual formulation for periodic solutions of systems undergoing frictional occurrences. 2023
preprint: https://hal.science/hal-04189699v1
    
------------------- Main variables -------------------
n - number of odd harmonics             mu - coefficient of friction
w - frequency                           Fn - normal force
m1, m2 - system masses                  rho - normalization parameter
k1, k2 - system stiffnesses
f1, f2 - system excitation forces

R, X_ - array of harmonic coefficients (Nonlinear force and positional quantities)
r, x_ - time series (Nonlinear force and positional quantities)
rdot, x_dot - time series (Nonlinear force and positional derivatives)

J1, J2 - conversion parameter for nonlinear force -> displacement (see below) (Ref: EQ. 22a, 22b)
F1, F2 - conversion parameter for nonlinear force -> displacement (see below) (Ref: EQ. 22a, 22b)

drdc, dxdc - piecewise jacobian coefficients (Ref: EQ. 28, 29)

Key: X1 = J1@R + F1
     X2 = J2@R + F2
     
------------------- Other variables ------------------
nval - (integer) harmonic coefficients array
num, dt, t, t_array - time parameters: number of steps, step size, time array (0,T), and array of harmonic timeseries (~IFFT procedure)

j_array, f_array - array of nonlinear force -> displacement conversion parameters (J1,J2 and F1,F2)
d_array - array of jacobian coefficients (drdc, dxdc)

GE - governig equation of motion f(R,x2dot, Fn, mu, rho) (Ref: EQ. 22c)
GE_J - jacobian of the the governing equation of motion f'(R,x2dot, Fn, mu, rho)

---------------------- Classes: ----------------------
solver - solver parameters for rootfinding procedure of governing equation (ref: EQ. 22c)
param - mechanical parameters of the 2DOF system 

Warning: viscous damping is not included!
"""
#%% Import relevant packages
import numpy as np
import time
from scipy.optimize import root

# Import relevant functions from the fun_Directory.py file
import fun_Directory # Custom functions used below
fun_FT = fun_Directory.fun_FT
fun_NonlinSystem = fun_Directory.fun_NonlinSystem
fun_Jac = fun_Directory.fun_Jac
fun_Plot = fun_Directory.fun_Plot

# Import class/structure from the class_Directory.py file
import class_Directory # Custom class for structures
SolverClass = class_Directory.SolverClass
ParamClass = class_Directory.ParamClass

#%% Main script
# System parameters
n = 100 # Number of harmonics
rho = 1 # parameter in friction equation
num = 1024 # Number of time-steps in one period
tols = 1e-6 # Nonlinear solver tolerance

# Parameters of the mechanical system
m1 = 1; m2 = 1
k1 = 1; k2 = 1
f1 = 0; f2 = 10
mu = 0.9; Fn = 8

# Forcing frequency and period
w = 0.4
T = 2*np.pi/np.array([w])

# Time array
t = np.linspace(0, T * (1 - 1 / 1024), 1024)
dt = t[1] - t[0] # time-step

# Initialize classes/structures/dictionnaries
Solver = SolverClass(None, None, None, n, None, num, dt, t, None, rho, w, Fn, mu, None, None, None)
Params = ParamClass(m1,m2,k1,k2,f1,f2)

# Construction of the required matrices and reduction to the unknowns of the problem (harmonics of the friction force)
fun_FT(Solver, Params)
# Initial guess R0 of the friction force harmonics
R0 = np.zeros(2*n)
# nonlinear solver
start_time = time.time()
result = root(fun_NonlinSystem, R0, args = Solver, method = 'hybr', jac=fun_Jac, tol = tols)
end_time = time.time()
elapsed_time = end_time - start_time

# Print and postprocess solver results
Solver.R = result.x # Result storage
print(f"Solver result: {result.success}"); 
print(f"Elapsed time: {elapsed_time} seconds")
fun_Plot(Solver)