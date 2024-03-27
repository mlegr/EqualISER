% Description: Implicit Euler for the forced-response for 2-dof oscillator
% with dry friction. Nonlinear friction is "approximated" via a
% prediction-correction scheme. At every time-step, sticking is assumed. If
% the resulting solution lies outside Coulomb's friction cone, sliding is
% enforced and the solution is updated.
% Author: Mathias Legrand
% Date: July 2, 2023

clear
clc

% mechanical parameters
run("../parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time domain solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1e-4;        % Time-step (h = 1e-4 is needed for convergence)
Nstep = 2e6;     % Total number of time steps. Sets the final time of the simulation
div = 1e3;       % Time interval for storage purposes. Storage every "div" time steps
Naf = Nstep/div; % Number of storage instances

% Initialisation of the solution array
sol = zeros(5,Naf); %(x1, v1, x2, v2, r)
sol_n = zeros(5,1);
sol_np = zeros(5,1);

% Sticking dynamic matrix: x2n+1 = x2n and v2n+1 = 0
% Unknowns: x1n+1, v1n+1, rn+1
Af=[1 -h 0;
   (k1+k2)*h m1+h*(d1+d2) 0;
   k2 d2 1];

% Sliding dynamic matrix: unknowns: x1n+1, v1n+1, x2n+1, v2n+1
As=[1 -h 0 0;
    (k1+k2)*h m1+h*(d1+d2) -k2*h -d2*h;
     0 0 1 -h;
    -h*k2 -h*d2 k2*h m2+h*d2];

% Forcing frequency and period of motion
w = omega;
T = 2*pi/w;
nT = floor(T/h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-domain loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time array
time = zeros(Naf,1);
t = 0;
tic
for n = 1 : Nstep
    t = t + h;

    % sol_n = sol(n) and sol_np = sol(n+1)
    x1 = sol_n(1,1);
    v1 = sol_n(2,1);
    x2 = sol_n(3,1);
    v2 = sol_n(4,1);

    % Assume sticking: x2n+1 = x2n and v2n+1 = 0.
    % Only three unknowns: x1, v1, r
    % Right hand-side dynamic vector
    bf = [x1;m1*v1+h*(f1*cos(w*t)+k2*x2);k2*x2-m2/h*v2-f2*cos(w*t)];
    new_state = Af\bf; % solve Af * new_state = bf
    rt = new_state(3,1); % get the friction force

    % Correction procedure
    % test if sticking is correct
    if abs(rt) < mu*N % sticking occurs
       sol_np(1,1) = new_state(1,1); % x1
       sol_np(2,1) = new_state(2,1); % v1
       sol_np(3,1) = sol_n(3,1);     % x2
       sol_np(4,1) = 0;              % v2
       sol_np(5,1) = rt;             % r
    else % sliding occurs. Four unknowns: x1, v1, x2, v2
       % right hand-side dynamic vector
       rt = mu*N*sign(rt); % should be: -mu*N*sign(v2) but v2 at the
       % on-going time step is unknown either (1) use a nonlinear solver on the friction law or (2)
       % assume that from sticking to sliding, the sliding force will be of
       % the sign of the sticking force at the previous time-step as done here
       bs = [x1; m1*v1+h*f1*cos(w*t); x2; m2*v2+h*(rt+f2*cos(w*t))];
       new_state = As\bs;
       sol_np(1:4,1) = new_state(1:4,1);
       sol_np(5,1) = rt;
    end

     % Storage
    counter = n/div;
    if (round(counter) - counter) == 0
        time(counter) = t;
        sol(:,counter) = sol_np(:,1);
        counter
    end 
    sol_n = sol_np;
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_PredCorr