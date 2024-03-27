% Description: HBM forced-response for 2-dof oscillator with Coulomb's friction on the
% second dof. Weighted-Residual formulation with true quadrature scheme
% This is much slower than the FFT version mostly because of the recursive
% quadrature scheme use to compute the integrals.
% Author: Mathias Legrand
% Date: June 10, 2023

clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mechanical parameters and matrices
run("../parameters")

% Number of harmonics (only odd harmonics). When Nphi is increased, it is
% beneficial to rely on a previous solution with fewer harmonics as initial
% guess of the nonlinear solver.
Nphi = 10;
SizeBlock = 2*Nphi; % cos + sin 

% external force (vector f) : cos 1 on first dof
force1 = zeros(SizeBlock,1);
force1(1,1) = f1; % f1*cos(wt) term on dof 1

% external force (vector f) : cos 1 on second dof
force2 = zeros(SizeBlock,1);
force2(1,1) = f2; % f2*cos(wt) term on dof 2 subject to friction

% Fourier/Frequency domain (structural) matrices
A11 = zeros(SizeBlock);
A22 = zeros(SizeBlock);

% Diagonal terms
% includes viscous damping
for n = 1:Nphi % cosine + sine: Only odd harmonics = 1, 3, 5, 7
    % dof 1
    A11(2*n-1,2*n-1) = k1 + k2 - (2*n-1)^2 * omega^2 * m1;
    A11(2*n-1,2*n) = (2*n-1) * omega * (d1+d2);
    A11(2*n,2*n) = k1 + k2 - (2*n-1)^2 * omega^2 * m1;
    A11(2*n,2*n-1) = -(2*n-1) * omega * (d1+d2);

    % dof 2
    A22(2*n-1,2*n-1) = k2 - (2*n-1)^2 * omega^2 * m2;
    A22(2*n-1,2*n) = (2*n-1) * omega * d2;
    A22(2*n,2*n) = k2 - (2*n-1)^2 * omega^2 * m2;
    A22(2*n,2*n-1) = -(2*n-1) * omega * d2;
end

% Off-diagonal terms
A12 = -k2*eye(SizeBlock,SizeBlock);
A21 = -k2*eye(SizeBlock,SizeBlock);
for n = 1:Nphi
    A21(2*n-1,2*n) = -(2*n-1) * omega * d2;
    A21(2*n,2*n-1) = (2*n-1) * omega * d2;
    A12(2*n-1,2*n) = -(2*n-1) * omega * d2; 
    A12(2*n,2*n-1) = (2*n-1) * omega * d2;
end

% Some temporary matrices and vectors
invD = inv(A22 - A21*(A11\A12)); % B21*(A1\B12) is equivalent to B21*inv(A1)*B12
force1on2 = -A21*(A11\force1); % force1 acting on x2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution initial guess
X0 = zeros(SizeBlock,1);

% Functions adot and phidot [works]
param=[mu,N,omega];

% Anonymous function for the system of nonlinear equations induced by
% friction. X = vector of harmonics of friction force r(t)
% With explicit Jacobian
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'TolFun',...
    1e-4,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
fun = @(X) fun_NonlinSystem_Jacobian(X,Nphi,param,invD,force1on2,force2);

% Without explicit Jacobian
%options = optimset('Algorithm','trust-region-dogleg','TolFun',1e-4,...
% 'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');%
%fun = @(X) fun_NonlinSystem(X,Nphi,param,invD,force1on2,force2);

tic
[r,fval,exitflag,output] = fsolve(fun,X0,options);
toc

% dof Frequency harmonics from friction force solution in Frequency domain
% (sol)
x2 = invD*(r + force2 + force1on2); %dof 2 in terms of sol
x1 = A11\(force1 - A12*x2); % dof 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_Fun