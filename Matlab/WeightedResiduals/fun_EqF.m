function F = fun_EqF(t,Coeff,Nphi,param,invD,f1on2,f2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transforms the velocity of second dof x2'(t) and friction force
% from frequency domain to time domain and computes the nonlinear friction
% force F in time domain.
% Is called by fun_NonlinSystem_Jacobian.m/fun_NonlinSystem.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t: time instant
% Coeff: array of friction force harmonics (in frequency domain)
% Nphi: number of (odd) harmonics in solution
% param: array of parameters
% invD: matrix to map friction force harmonics to dof 2 harmonics
% f1on2: effect of force 1 on dof 2 in frequency domain
% f2: force 2 on dof 2 in frequency domain

omega=param(1,3);
mu=param(1,1);
N=param(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectorized cos and sin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_p = omega.*(1:2:2*Nphi-1).';
cot = cos(omega_p.*t);
sot = sin(omega_p.*t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friction force r(t) in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_t = sum(Coeff(1:2:end-1,1).*cot + Coeff(2:2:end,1).*sot,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocity of second dof x2'(t) (retrieved from Coeff) in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = invD*(f2 + f1on2 + Coeff); % b in terms of c
xdot_t = sum(omega_p.*(b(2:2:end,1).*cot - b(1:2:end-1,1).*sot),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear function (in time domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Stadler Eq 2.6)
rho = 1;
F = xdot_t+min(0,rho*(mu*N+r_t)-xdot_t)+max(0,rho*(-mu*N+r_t)-xdot_t); %Eq (17c)

% C. Pierre 1985. Does not work when sticking occurs/ no manual Jacobian
% Output = r_t+sign(xdot_t); 

% Berthillier 1998 ASME/no manual Jacobian
%delta = 0.125;
%Output = mu*N*xdot_t/delta+r_t+min(mu*N*(1-xdot_t/delta),0)+max(-mu*N*(1+xdot_t/delta),0); 
