% solution plot over one period.
% Is called by IntFreq.m

% function plots [works]
t = linspace(0,2*pi/omega,2000);

fric = 0*ones(size(t))/2;

x2_t = 0*ones(size(t))/2;
v2_t = zeros(size(t));

x1_t = 0*ones(size(t))/2;
v1_t = zeros(size(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friction force, displacements and velocities in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:Nphi
    omega_p = (2*p-1)*omega;
    fric = fric + r(2*p-1,1).*cos(omega_p.*t)+r(2*p,1).*sin(omega_p.*t);

    % dof 2
    x2_t = x2_t + x2(2*p-1,1).*cos(omega_p.*t)+x2(2*p,1).*sin(omega_p.*t);
    v2_t = v2_t + omega_p*(x2(2*p,1).*cos(omega_p.*t)-x2(2*p-1,1).*sin(omega_p.*t));
    
    % dof 1
    x1_t = x1_t + x1(2*p-1,1).*cos(omega_p.*t) + x1(2*p,1).*sin(omega_p.*t);
    v1_t = v1_t + omega_p*(x1(2*p,1).*cos(omega_p.*t) - x1(2*p-1,1).*sin(omega_p.*t));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% main title
sgtitle(['N_{\phi}=', num2str(Nphi), ', \mu = ', num2str(mu), ', N = ', num2str(N), ', \omega = ', num2str(omega)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,1)
plot(t,x2_t,'k',t,x1_t,'b')
ylabel('displacement')
legend('x2','x1','Location','southeast')
xlabel('time')
axis tight

subplot(2,3,2)
plot(t,v2_t,'k',t,v1_t,'b')
ylabel('velocity')
legend('v2','v1','Location','southeast')
xlabel('time')
axis tight

subplot(2,3,3)
plot(t,fric,'r')
xlabel('time')
ylabel('r')
axis tight

subplot(2,3,4)
plot(v2_t,fric,'k')
ylabel('r(t)')
xlabel('v2(t)')
axis tight

subplot(2,3,5)
plot(x2_t,fric,'k')
ylabel('r(t)')
xlabel('x2(t)')
axis tight

subplot(2,3,6)
plot(x1_t,fric,'k')
ylabel('r(t)')
xlabel('x1(t)')
axis tight