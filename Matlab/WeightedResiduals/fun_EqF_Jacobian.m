function J = fun_EqF_Jacobian(t,Coeff,Nphi,param,invD,f1on2,f2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates the Jacobian matrix J of the nonlinear function F computed in
% fun_EqF.m.
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
rho = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectorized cos and sin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_p = omega.*(1:2:2*Nphi-1).';
cot = cos(omega_p.*t);
sit = sin(omega_p.*t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Friction force r(t) in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_t = sum(Coeff(1:2:end-1,1).*cot + Coeff(2:2:end,1).*sit,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity of second dof x2'(t) (retrieved from Coeff) in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=invD*(f2+f1on2+Coeff); % b in terms of c
xdot_t = sum(omega_p.*(b(2:2:end,1).*cot - b(1:2:end-1,1).*sit),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian (in time domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test on which term is active in the nonlinear function 
% xdot_t+min(0,rho*(mu*N+r_t)-xdot_t)+max(0,rho*(-mu*N+r_t)-xdot_t)
% See piecewise definition of the nonlinear function in the paper

Output = zeros(2*Nphi,1); % Eq (29)
if rho*(r_t+mu*N)-xdot_t > 0 && rho*(r_t-mu*N)-xdot_t < 0
    Output(1:2:end-1,1) = - omega_p.*sit;
    Output(2:2:end,1) = omega_p.*cot;
    J = Output'*invD;
else
    Output(1:2:end-1,1) = rho * cot;
    Output(2:2:end,1) = rho * sit;
    J = Output';
end