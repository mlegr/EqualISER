function [F, J] = fun_NonlinSystem(Coeff,Nphi,param,invD,f1on2,f2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates the system of nonlinear equations F = 0 to be solved for
% and attendant Jacobian matrix J.
% Calls fun_R.m and fun_xdot.m
% Is called by main.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeff: array of friction force harmonics (in frequency domain)
% Nphi: number of (odd) harmonics in solution
% param: array of parameters (see main.m)
% invD: matrix to map friction force harmonics to dof 2 harmonics  (see main.m)
% f1on2: effect of force 1 on dof 2 in frequency domain (see main.m)
% f2: force 2 on dof 2 in frequency domain (see main.m)

% Retrieve forcing frequency omega, friction coefficient mu and normal
% force N
omega = param(1,3);
mu = param(1,1);
N = param(1,2);
rho = 1; % rho parameter in friction equality

% period
T=2*pi/omega;
omega_p = omega.*(1:2:2*Nphi-1)';

nt = 1024; % number of time steps in the FFT (can be increased in powers of 2)
time = linspace(0,T*(1-1/nt),nt);
r_t = fun_R(Coeff,nt); % friction force in time domain
xdot_t = fun_xdot(Coeff,Nphi,param,invD,f1on2,f2,nt); %x2' in time domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear friction equation ψr in time domain
% Stadler (2004). SIAM Journal on Optimization 15(1):39–62, DOI 10.1137/S1052623403420833,
EqFval = xdot_t + min(0,rho*(mu*N + r_t) - xdot_t) + max(0,rho*(-mu*N + r_t) - xdot_t); % (Stadler Eq 2.6)

% FFT on ψr: equivalent to integral or weighted residuals with Fourier functions
vec = 2*fft(EqFval)/nt;
F = zeros(2*Nphi,1);
F(1:2:end-1,1) = real(vec(2:2:2*Nphi));
F(2:2:end,1) = -imag(vec(2:2:2*Nphi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corresponding Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialiation of ∂ψr/∂r in time domain
dpsi_dr = zeros(2*Nphi,nt);

% Test on the portions of the nonlinear function to build ∂ψr/∂r
idx = r_t + mu*N - xdot_t > 0 & r_t - mu*N - xdot_t < 0;
dpsi_dr(1:2:end-1,idx) = - omega_p.*sin(omega_p.*time(idx));
dpsi_dr(2:2:end,idx) = omega_p.*cos(omega_p.*time(idx));
dpsi_dr = invD'*dpsi_dr;
dpsi_dr(1:2:end-1,~idx) = rho * cos(omega_p.*time(~idx));
dpsi_dr(2:2:end,~idx) = rho * sin(omega_p.*time(~idx));

% FFT on ∂ψr/∂r: equivalent to integral or weighted residuals with Fourier functions
mat_J = 2*fft(dpsi_dr,nt,2)/nt;
J = zeros(2*Nphi,2*Nphi);
J(1:2:end-1,:) = real(mat_J(:,2:2:2*Nphi))';
J(2:2:end,:) = -imag(mat_J(:,2:2:2*Nphi))';