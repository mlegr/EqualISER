function [F] = fun_NonlinSystem(Coeff,Nphi,param,invD,f1on2,f2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates the system of nonlinear equations F to be solved for.
% Calls fun_EqF.m.
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

% period
omega = param(1,3);
T = 2*pi/omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tolerances of the quadrature schemes
optsF = {'AbsTol',1.e-4,'RelTol',0};

fun = @(t) [cos((2*(1:Nphi)-1).*omega.*t),sin((2*(1:Nphi)-1).*omega.*t)]'.*fun_EqF(t,Coeff,Nphi,param,invD,f1on2,f2);
F = integral(fun,0,T,optsF{:},'ArrayValued',true);