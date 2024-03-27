function F = fun_NonSmoothDynCoulomb(var,param)

% Sets the system of nonlinear equations F = 0 to be soved for using Implicit Euler
% on the governing equations

t=param(1);
mu=param(2);
N=param(3);
h=param(4);
m1=param(5);
m2=param(6);
k1=param(7);
k2=param(8);
f1=param(9);
f2=param(10);
x1p=param(11);
v1p=param(12);
x2p=param(13);
v2p=param(14);

w=param(15);

x1=var(1);
v1=var(2);
x2=var(3);
v2=var(4);
r=var(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsmooth nonlinear equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear equations
F(1)=x1-h*v1-x1p;
F(2)=m1*(v1-v1p)-h*f1*cos(w*t)+h*((k1+k2)*x1-k2*x2);
F(3)=x2-h*v2-x2p;
F(4)=m2*(v2-v2p)-h*f2*cos(w*t)+h*k2*(x2-x1)-h*r;

% Nonlinear friction equation
rho = 1;
F(5) = v2 + min(0,rho*(r + mu*N) - v2) + max(0,rho*(r - mu*N) - v2);