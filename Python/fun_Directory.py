"""
Title: FD Solver for ML-2DOF Oscillator with friction
Author: Sebastian McDonald
Date: Feb 29, 2024    
    
Description: Function directory for FD solver for pure friction nonlinearity of 2DOF oscillator
Reference: Mathias Legrand, Christophe Pierre. A compact, equality-based weighted residual
formulation for periodic solutions of systems undergoing frictional occurrences. 2023. ⟨hal-04189699⟩
Viscous damping is ignored.

----------------- Function definition -----------------
fun_FT - initialize the conversion parameters and jacobian coefficients (J1, J2, F1, F2, dxdc, and drdc) 
fun_NonlinSystem - governing equation of motion: f(R,x2dot, Fn, mu, rho) (Ref: EQ. 22c)
fun_x2dot - transformaton function: R -> X2 -> x2dot
fun_Jac - jacobian of fun_NonlinSystem
fun_Plot - postprocessing function

"""
#%% Import relevant packages
import numpy as np
from scipy.fft import fft
import matplotlib.pyplot as plt

#%% fun_FT
def fun_FT(Solver, Param):
    #  Unpack object classes
    n = Solver.n
    w = Solver.w
    t = Solver.t
    m1 = Param.m1
    m2 = Param.m2
    k1 = Param.k1
    k2 = Param.k2
    f1 = Param.f1
    f2 = Param.f2

    # Initialize and define conversion parameters (J1, J2)
    nval = np.array([2*i + 1 for i in range(n)])
    nsquare = nval**2
    J1 = np.zeros((2*n, 2*n))
    J2 = np.zeros((2*n, 2*n))
    F1 = np.zeros((2*n,1))
    F2 = np.zeros((2*n,1))
    F11 = 0
    F21 = 0
    
    denom = (k1*k2 - nsquare*k2*m1*(w**2) - nsquare*k1*m2*(w**2) - nsquare*k2*m2*(w**2) + (nsquare**2)*(m1*m2)*(w**4))
    
    jth1 = np.array([[np.array([k2]), np.array([0])], [np.array([0]), np.array([k2])]])/denom  #stored in 3rd axis
    jth2 = np.array([[k1 + k2 - nsquare*m1*w**2, np.zeros(len(nsquare))], [np.zeros(len(nsquare)), k1 + k2 - nsquare*m1*w**2]])/denom
    
    for i in range(len(nval)):
        J1[nval[i]-1:(nval+1)[i], nval[i]-1:(nval+1)[i]] = jth1[:,:,i]
        J2[nval[i]-1:(nval+1)[i], nval[i]-1:(nval+1)[i]] = jth2[:,:,i]
    
    denomf = denom[0:1]
    F11 = np.array([[f1*k2 + f2*k2  - f1*m2*w**2],[0]])/denomf[:, np.newaxis]
    F21 = np.array([[f2*k1 + f1*k2 + f2*k2 - f2*m1*w**2],[0]])/denomf[:, np.newaxis]
    F1[0:2] = F11
    F2[0:2] = F21
    
    # Evaluate the harmonic time series arrays (equivalent to iFFT)
    w = np.array([w]); i = nval
    iw = np.einsum('i, k-> ik',i,w)
    iwt = np.einsum('i, jk, k-> jik',i,t,w)
    array1 = np.zeros((len(t), 2*n, np.size(w)))
    array2 = np.zeros((len(t), 2*n, np.size(w)))
    array3 = np.zeros((len(t), 2*n, np.size(w)))
    
    sin_values1 = np.sin(iwt);
    cos_values1 = np.cos(iwt);
    array1[:, 0::2, :] = cos_values1
    array1[:, 1::2, :] = sin_values1
    
    sin_values2 = -(iw)*np.sin(iwt)
    cos_values2 = (iw)*np.cos(iwt)
    array2[:, 0::2, :] = sin_values2
    array2[:, 1::2, :] = cos_values2
    
    sin_values3 = -iw**2*np.sin(iwt)
    cos_values3 = -iw**2*np.cos(iwt)
    array3[:, 0::2, :] = cos_values3
    array3[:, 1::2, :] = sin_values3
        
    # Evaluate the piecewise jacobian coefficients
    drdc = np.squeeze(array1).T
    dxdc = np.einsum('ijk,kj->ji', array2, J2.T)
    
    # Write variables to object class Solver
    Solver.nval = nval
    Solver.j_array = [J1, J2]
    Solver.f_array = [F1, F2]
    Solver.t_array = [array1, array2, array3]
    Solver.d_array = [drdc, dxdc]

    return

#%% fun_NonlinSystem
def fun_NonlinSystem(R0, Solver):
    # Unpack object classes
    R0 = np.atleast_2d(R0)
    Solver.R = R0
    w = Solver.w
    n = Solver.n
    rho = Solver.rho
    Fn = Solver.Fn
    mu = Solver.mu
    array1 = Solver.t_array[0]
    num = Solver.num
    
    # Evaluate r and x2dot
    r = np.einsum('ijk,kj->ki', array1, R0)
    Solver.r = r
    x2dot = fun_x2dot(Solver) 
    Solver.x2dot = x2dot
    
    # Evaluate the governing nonsmooth equation f(R,x2dot, Fn, mu, rho) 
    GE = x2dot + np.minimum(0, rho*(r + mu*Fn) - x2dot) + np.maximum(0, rho*(r - mu*Fn) - x2dot)
    
    # Compute R via FFT (nonlinear/nonsmooth function in frequency domain)
    Rvec = fft(GE)/(num/(2*np.pi/w))
    R = np.zeros((2*n,1))
    R[0::2, 0] = np.real(Rvec[0, 1:2*n + 1:2])
    R[1::2, 0] = -np.imag(Rvec[0, 1:2*n + 1:2])
    
    return np.squeeze(R)

#%% fun_x2dot
def fun_x2dot(Solver):
    # Unpack object classes
    R = Solver.R
    J2 = Solver.j_array[1]
    F2 = Solver.f_array[1]
    array2 = Solver.t_array[1]
    
    # Transform R -> X2
    X2 = np.einsum('ii,ij->ij', J2, R.T) + F2
    
    # Transform X2 -> X2dot (Equivalent to the iFFT procedure)
    x2dot = np.einsum('ijk,kj->ki', array2, X2.T)
    
    return np.squeeze(x2dot)

#%% fun_Jac (Computation of the Jacobian)
def fun_Jac(R0, Solver):
    # Unpack object classes
    w = Solver.w
    num = Solver.num
    n = Solver.n
    rho = Solver.rho
    Fn = Solver.Fn
    mu = Solver.mu
    drdc = Solver.d_array[0]
    dxdc = Solver.d_array[1]
    
    # Note: rewriting arrays in jacobian seems to be moderately faster than parsing them from SolverClass
    #r = Solver.r; x2dot = Solver.x2dot; 
    R0 = np.atleast_2d(R0)
    array1 = Solver.t_array[0]
    r = np.einsum('ijk,kj->ki', array1, R0)
    x2dot = fun_x2dot(Solver)
    
    # Compute the indices (idx) associated with the piecewise derivative (GE_jac)
    idx1 = np.where(((rho*(r + mu*Fn) - (x2dot)) > 0) & ((rho*(r - mu*Fn) - (x2dot)) > 0))[1]
    idx2 = np.where(((rho*(r + mu*Fn) - (x2dot)) > 0) & ((rho*(r - mu*Fn) - (x2dot)) < 0))[1]
    idx3 = np.where(((rho*(r + mu*Fn) - (x2dot)) < 0)  & ((rho*(r - mu*Fn) - (x2dot)) < 0))[1]
  
    GE_jac = np.zeros((2*n, num))
    
    GE_jac[:,idx1] = rho*(drdc[:,idx1])
    GE_jac[:,idx2] = 1*(dxdc[:,idx2])
    GE_jac[:,idx3] = rho*(drdc[:,idx3])
    
    # Compute R_J via FFT
    Rmat_j = fft(GE_jac)/(num/(2*np.pi/w))
    R_J = np.zeros((2*n,2*n))
    R_J[0::2, :] = np.real(Rmat_j[:, 1:2*n + 1:2]).T
    R_J[1::2, :] = -np.imag(Rmat_j[:, 1:2*n + 1:2]).T
  
    return R_J

#%% fun_Plot
def fun_Plot(Solver):
    # Unpack object class solver
    R = Solver.R
    R = np.atleast_2d(R)
    w = Solver.w
    w = np.array([w])
    t = Solver.t
    n = Solver.n
    J1 = Solver.j_array[0]
    J2 = Solver.j_array[1]
    F1 = Solver.f_array[0]; F2 = Solver.f_array[1]
    array1 = Solver.t_array[0]
    # array2 = Solver.t_array[1] (see below for the velocity)

    # Evaluate quantities of interest via appropriate transformation procedures
    X1 = J1@R.T + F1
    X2 = J2@R.T + F2
    x1 = np.einsum('ijk,kj->ki', array1, X1.T)
    x2 = np.einsum('ijk,kj->ki', array1, X2.T)
    r = np.einsum('ijk,kj->ki', array1, R)
    # x2dot = np.einsum('ijk,kj->ki', array2, X2.T) # velocity in case it is to be plotted

    # Postprocess and plot results - currently only showing R (the friction force)
    xplot = np.array([x1,x2])
    tplot = np.array([t,t])
    xplot
    tplot
    plt.plot(t, r.T)
    plt.title("Nonlinear force plot (n = {})".format(n))
    plt.xlabel("Time (s)")
    plt.ylabel("Newtons (N)")    
    return