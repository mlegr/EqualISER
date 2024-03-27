function [Output] = fun_NonlinSystem_NumQuad(Coeff,Nphi,param,DiaginvD,f1,f2,ScaledInt,dif1,wgt1)

% Replaces the built-in quadrature scheme of Matlab ('integral' command) by a 
% custom numerical scheme. 
% This is just an attempt.

omega=param(1,3);
Output=zeros(2*Nphi,1);
for p = 1 : Nphi

    fun1 = @(t) cos((2*p-1)*omega*t);
    fun2 = @(t) sin((2*p-1)*omega*t);
    fun3 = @(t) fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2);

    % Evaluate f at scaled intervals, need to map each interval to [-1 1].
    an1 = fun1(ScaledInt);      % Function evaluations.
    an2 = fun2(ScaledInt);      % Function evaluations.
    an3 = arrayfun(fun3,ScaledInt);    % Function evaluations.
   % an4 = arrayfun(fun3,ScaledInt);   % Function evaluations.
    new1 = dif1*(an1.*an3)'*wgt1;      % Multiply by the weights and intervals.
    new2 = dif1*(an2.*an3)'*wgt1;                % Multiply by the weights and intervals.
    Output(2*p-1,1) = sum(new1(:));
    Output(2*p,1) = sum(new2(:));

    % Galerkin on EqF
    %Output(2*p-1,1) = integral(@(t) cos((2*p-1)*omega*t).*mem_fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
    %Output(2*p,1) = integral(@(t) sin((2*p-1)*omega*t).*mem_fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
    %Output(2*p-1,1) = integral(@(t) cos((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
    %Output(2*p,1) = integral(@(t) sin((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
    %Output(2*p-1,1) = gaussquad(@(t) cos((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T);
    %Output(2*p,1) = gaussquad(@(t) sin((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T);
end