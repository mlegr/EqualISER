% Description: sets the mechanical system of interest with other various parameters
% Is called by the main scripts in other folders

% Mechanical parameters
Ndof = 2; % Number of dof
m1 = 1; k1 = 1; % masses and stiffnesses
m2 = 1; k2 = 1;
d1 = 0.02; d2 = 0.02;% Visvcous damping

% Friction 
mu = 0.9; N = 9;

% External forces
f1 = 20; % force magnitude on dof 1
f2 = 0;  % force magnitude on dof 2

% Number of harmonics (in the Frequency Domain solvers)
Nphi = 10;

% Forcing frequency in F*cos(ometa*t)
omega = .4;

% Structural Matrices
M = [m1,0;0,m2];
K = [k1+k2,-k2;-k2,k2];
D = [d1+d2,-d2;-d2,d2];