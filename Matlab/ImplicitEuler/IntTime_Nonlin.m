% Description: Implicit Euler for the forced-response for 2-dof oscillator
% with dry friction. Nonlinear friction is truly solved and not approximated. 
% Much slower than IntTimePredCorr.m. In the current version, the viscous
% damping is ignored, that is d1 = d2 = 0 should be set in parameters.m for
% comparison purposes.
% Author: Mathias Legrand
% Date: July 2, 2023

clear
clc

% mechanical parameters
run("../parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time domain solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1e-4;        % Time-step  (h = 1e-4 is needed for convergence)
Nstep = 2e6;     % Total number of time steps. Sets the final time of the simulation
div = 1e3;       % Time interval for storage purposes. Storage every "div" time steps
Naf = Nstep/div; % Number of storage instances

% Initialisation of the solution array
sol = zeros(5,Naf); %(x1, v1, x2, v2, r)

% Initialization
x1 = 0;
x2 = 0;
v1 = 0;
v2 = 0;
r = 0;
X0 = zeros(5,1);

% array of parameters (to be used as input below)
param = zeros(15,1);
param(2) = mu;
param(3) = N;
param(4) = h;
param(5) = m1;
param(6) = m2;
param(7) = k1;
param(8) = k2;
param(9) = f1;
param(10) = f2;
param(15) = omega;

% Forcing frequency and period of motion
w = omega;
T = 2*pi/w;
nT = floor(T/h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-domain loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = zeros(Naf,1);
t = 0;
tic
for n = 1 : Nstep
    t = t + h;

    % Nonlinear equation set-up
    param(1) = t;
    param(11) = x1;
    param(12) = v1;
    param(13) = x2;
    param(14) = v2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nonlinear solver 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Without explicit Jacobian: comment/uncomment as needed
    options = optimset('TolFun',1e-6,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off');
    fun = @(X) fun_NonSmoothDynCoulomb(X,param);

    % With explicit Jacobian: comment/uncomment as needed.
    % Much slower than above: not clear why.
    % fun = @(X) fun_NonSmoothDynCoulomb_Jacobian(X,param); %Works! Init
    % options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'TolFun',1e-3,'MaxFunEvals',...
	%	1e5,'Maxiter',1e5,'Display','none');
    
    [solp] = fsolve(fun,X0,options);
    X0 = solp; % Initial guess for next time step

    % Update for next iteration
    x1 = solp(1);
    v1 = solp(2);
    x2 = solp(3);
    v2 = solp(4);

     % Storage
    counter = n/div;
    if (round(counter) - counter) == 0
        time(counter) = t;
        sol(:,counter) = solp(:,1);
        counter
    end 
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_Nonlin