function [Output] = fun_NonlinSystem_NumQuad_New(Coeff,Nphi,param,invD,f1,f2)

% Replaces the built-in quadrature scheme of Matlab ('integral' command) by a 
% custom numerical scheme. 
% This is just an attempt.

omega=param(1,3);
T = 2*pi/omega;

%opts = {'AbsTol',1.e-6,'RelTol',0}; %
%opts = {'AbsTol',1.e-6,'RelTol',1.e-6};

% Initially use a 40 point Gauss-Legendre quadrature in 15 subintervals.
gp = 20; %40  
ints = 5;  %15
[abs1, wgt1] = Gauss(gp);                               % Get Gauss points.
bb(1) = 0;                                              % Get subintervals.
bb(2:ints) = (2:ints)*T/ints;
dif = diff(bb)/2;
ScaledInt = (abs1+1)*dif+repmat(bb(1:end-1),gp,1);
% harmonics in phi (cos and sin)

% hallo = @(t) fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2);
% mem_fun_EqF = memoize(@fun_EqF);
Output=zeros(2*Nphi,1);

fun3 = @(t) fun_EqF(t,Coeff,Nphi,param,invD,f1,f2);
an3 = arrayfun(fun3,ScaledInt);      % Function evaluations.

for p = 1 : Nphi
	fun1 = @(t) cos((2*p-1)*omega*t);
	fun2 = @(t) sin((2*p-1)*omega*t);

	% Evaluate f at scaled intervals, need to map each interval to [-1 1].
	an1 = fun1(ScaledInt);      % Function evaluations.
	an2 = fun2(ScaledInt);      % Function evaluations.
	%an4 = arrayfun(fun3,ScaledInt);      % Function evaluations.
	new1 = dif*(an1.*an3)'*wgt1;                % Multiply by the weights and intervals.
	new2 = dif*(an2.*an3)'*wgt1;                % Multiply by the weights and intervals.
	Output(2*p-1,1) = sum(new1(:));
	Output(2*p,1) = sum(new2(:));

	% Galerkin on EqF
	% Output(2*p-1,1) = integral(@(t) cos((2*p-1)*omega*t).*mem_fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
	% Output(2*p,1) = integral(@(t) sin((2*p-1)*omega*t).*mem_fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
	% Output(2*p-1,1) = integral(@(t) cos((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
	% Output(2*p,1) = integral(@(t) sin((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T, opts{:});
	% Output(2*p-1,1) = gaussquad(@(t) cos((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T);
	% Output(2*p,1) = gaussquad(@(t) sin((2*p-1)*omega*t).*fun_EqF(t,Coeff,Nphi,param,DiaginvD,f1,f2), 0, T);
end
end


function [x, w] = Gauss(n)
    % Generates the abscissa and weights for a Gauss-Legendre quadrature.
    % Reference:  Numerical Recipes in Fortran 77, Cornell press.
    x = zeros(n,1);                                           % Preallocations.
    w = x;
    m = (n+1)/2;
    for ii=1:m
        z = cos(pi*(ii-.25)/(n+.5));                        % Initial estimate.
        z1 = z+1;
        while abs(z-z1)>eps
            p1 = 1;
            p2 = 0;
            for jj = 1:n
                p3 = p2;
                p2 = p1;
                p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj;       % The Legendre polynomial.
            end
            pp = n*(z*p1-p2)/(z^2-1);                        % The L.P. derivative.
            z1 = z;
            z = z1-p1/pp;
        end
        x(ii) = -z;                                   % Build up the abscissas.
        x(n+1-ii) = z;
        w(ii) = 2/((1-z^2)*(pp^2));                     % Build up the weights.
        w(n+1-ii) = w(ii);
    end
end

