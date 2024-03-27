"""
Title: FD Solver for ML-2DOF Oscillator with friction
Author: Sebastian McDonald
Date: Feb 29, 2024    
    
Description: Class directory creating a structure-like data array in Python
"""

#%% Pythonic structures storing important parameters to be called within the various functions
# A dictionnary could have done the job as well
class SolverClass:
    def __init__(self, R, r, x2dot, n, nval, num, dt, t, t_array, rho, w, Fn, mu, j_array, f_array, d_array):
        self.R = R; 
        self.r = r; 
        self.x2dot = x2dot;
        
        self.n = n; 
        self.nval = nval; 
        
        self.num = num; 
        self.dt = t; 
        self.t = t; 
        self.t_array = t_array;

        self.rho = rho; 
        self.w = w;
        
        self.Fn = Fn; 
        self.mu = mu;
        
        self.j_array = j_array; 
        self.f_array = f_array
        
        self.d_array = d_array;

class ParamClass:
    def __init__(self, m1, m2, k1, k2, f1, f2):
        self.m1 = m1; self.m2 = m2;
        self.k1 = k1; self.k2 = k2;
        self.f1 = f1; self.f2 = f2;