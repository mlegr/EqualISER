function Output = fun_R(Coeff,nt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Friction force from frequency domain to time domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeff: array of friction force harmonics (in frequency domain)
% nt: number of time-steps in a period of motion


lCoeff = length(Coeff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Friction force in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Harm = zeros(nt,1);
Harm(2:2:lCoeff) = Coeff(1:2:lCoeff-1)-1i*Coeff(2:2:lCoeff);% start
Harm(nt:-2:nt-lCoeff+1) = Coeff(1:2:lCoeff-1)+1i*Coeff(2:2:lCoeff);% end
Output = ifft(Harm,nt)'*nt/2;