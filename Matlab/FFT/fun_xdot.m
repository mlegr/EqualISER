function Output = fun_xdot(Coeff,Nphi,param,invD,f1on2,f2,nt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trasnforms the velocity of second dof x2'(t) from frequency domain to
% time domain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeff: array of friction force harmonics (in frequency domain)
% Nphi: number of (odd) harmonics in solution
% param: array of parameters (see main.m)
% invD: matrix to map friction force harmonics to dof 2 harmonics  (see main.m)
% f1on2: effect of force 1 on dof 2 in frequency domain (see main.m)
% f2: force 2 on dof 2 in frequency domain (see main.m)
% nt: number of time-steps in a period of motion


% Retrieve forcing frequency
omega=param(1,3);
omega_p = omega.*(1:2:2*Nphi-1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocity of second dof x2'(t) in frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=invD*(f2+f1on2+Coeff); % b in terms of c
lCoeff = length(Coeff);
bCoeff = zeros(length(Coeff),1);
bCoeff(2:2:lCoeff) = -omega_p.*b(1:2:end);%sin terms
bCoeff(1:2:lCoeff-1) = omega_p.*b(2:2:end);%cos terms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocity of second dof x2'(t) in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Harm = zeros(nt,1);
Harm(2:2:lCoeff) = bCoeff(1:2:lCoeff-1)-1i*bCoeff(2:2:lCoeff);% start
Harm(nt:-2:nt-lCoeff+1) = bCoeff(1:2:lCoeff-1)+1i*bCoeff(2:2:lCoeff);% end
Output = ifft(Harm,nt)'*nt/2;