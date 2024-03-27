function [F] = fun_NonlinSystem_Vect(Coeff,Nphi,param,invD,f1on2,f2)

omega=param(1,3);

% period
T=2*pi/omega;

opts = {'AbsTol',1.e-4,'RelTol',0}; %
%opts = {'AbsTol',1.e-6,'RelTol',1.e-6};

% harmonics in phi (cos and sin)
%fun = @(t) [cos((2*(1:Nphi)-1).*omega.*t),sin((2*(1:Nphi)-1).*omega.*t)].*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2);
F = integral(@(t) [cos((2*(1:Nphi)-1).*omega.*t),sin((2*(1:Nphi)-1).*omega.*t)].*fun_EqF(t,Coeff,Nphi,param,invD,f1on2,f2),0,T,'ArrayValued',true, opts{:});

%Output=zeros(SizeBlock,1);
%for p = 1 : Nphi
    % Galerkin on EqF
    %Output(2*p-1,1) = integral(@(t) cos((2*p-1)*omega*t).*mem_fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
    %Output(2*p,1) = integral(@(t) sin((2*p-1)*omega*t).*mem_fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
 %   Output(2*p-1,1) = integral(@(t) cos((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
 %   Output(2*p,1) = integral(@(t) sin((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
    %Output(2*p-1,1) = gaussquad(@(t) cos((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T);
    %Output(2*p,1) = gaussquad(@(t) sin((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T);
%end