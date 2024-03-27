function [F,J] = fun_NonSmoothDynCoulomb_Jacobian(var,param)

% Sets the system of nonlinear equations F = 0 to be soved for using Implicit Euler
% on the governing equations, along with the Jacobian matrix J

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

x1=var(1); % displacement dof 1
v1=var(2); % velocity dof 1
x2=var(3); % displacement dof 2
v2=var(4); % velocity dof 2
r=var(5);  % friction force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsmooth nonlinear equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear equations
F(1) = x1-h*v1-x1p;
F(2) = m1*(v1-v1p)-h*f1*cos(w*t)+h*((k1+k2)*x1-k2*x2);
F(3) = x2-h*v2-x2p;
F(4) = m2*(v2-v2p)-h*f2*cos(w*t)+h*k2*(x2-x1)-h*r;

% nonlinear friction equation equality style
rho = 1;
F(5) = v2+min(0,rho*(r+mu*N)-v2)+max(0,rho*(r-mu*N)-v2); %Eq (17c)

% Corresponding Jacobian
H1 = heaviside(rho*(mu*N+r)-v2);
H2 = heaviside(rho*(r-mu*N)-v2);
J = zeros(5,5);
J(1,1) = 1; J(1,2)=-h;
J(2,1) = (k1+k2)*h; J(2,2) = m1; J(2,3) = -k2*h;
J(3,3) = 1; J(3,4) = -h;
J(4,1) = -h*k2; J(4,3) = h*k2; J(4,4) = m2; J(4,5) = -h;
J(5,4) = H1-H2; J(5,5) = rho*(1-H1+H2);