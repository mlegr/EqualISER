% Description: HBM forced-response for 2-dof oscillator with Coulomb's friction on the
% second dof. FFT is used in the quadrature to replace the classical
% and slower quadrature schemes. Similar to main.m but adapted to create a
% Frequency Response plot. Called by Freq_Plot.m
% Author: Mathias Legrand
% Date: March 22, 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of harmonics (only odd harmonics participate in solution)
SizeBlock = 2*Nphi; % cos + sin

% External force (vector f)
force1 = zeros(SizeBlock,1);
force1(1,1) = f1; % f1*cos(wt) term on dof 1

% external force (vector f)
force2 = zeros(SizeBlock,1);
force2(1,1) = f2; % f2*cos(wt) term on dof 2 (subject to friction)

% Param array as input below
param=[mu,N,omega];

% Fourier/Frequency domain (structural) matrices
A11 = zeros(SizeBlock);
A22 = zeros(SizeBlock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamics matrices in frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Temporary matrices and vectors
invD = inv(A22-A21*(A11\A12)); % B21*(A1\B12) is equivalent to B21*inv(A1)*B12
force1on2 = -A21*(A11\force1); % force1 acting on x2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solution initial guess
if kk == 1
    X0 = zeros(SizeBlock,1);
% Uses the previous number of harmonics in initial guess
else
    ns=length(r);
    X0(1:ns,1) = r;
end

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'TolFun',1e-6,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','none');
options.Algorithm = 'trust-region-dogleg';
fun = @(X) fun_NonlinSystem(X,Nphi,param,invD,force1on2,force2);

tic
[r,fval,exitflag,output] = fsolve(fun,X0,options);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
x2 = invD*(r+force2+force1on2); %dof 2 in terms of sol
x1 = A11\(force1-A12*x2); % dof 1
%Plot_Fun