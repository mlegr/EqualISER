% solution plot over one period. Run main.m file before
% close all; clc;
nt = 1024; % number of time-steps in the period of the motion/forcing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots: initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = linspace(0,T,nt);
fric = 0*ones(size(t)); % friction force

x2_t = 0*ones(size(t)); % dof 2 displacement
v2_t = zeros(size(t));    % dof 2 velocity

x1_t = 0*ones(size(t)); % dof 1 displacement
v1_t = zeros(size(t));    % dof 1 velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Friction force, displacement, velocity in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:Nphi
    omega_p = (2*p-1)*omega;

    % friction force
    fric = fric + r(2*p-1,1).*cos(omega_p.*t) + r(2*p,1).*sin(omega_p.*t);

    % dof 2
    x2_t = x2_t + x2(2*p-1,1).*cos(omega_p.*t) + x2(2*p,1).*sin(omega_p.*t);
    v2_t = v2_t + omega_p*(x2(2*p,1).*cos(omega_p.*t) - x2(2*p-1,1).*sin(omega_p.*t));
    
    % dof 1
    x1_t = x1_t + x1(2*p-1,1).*cos(omega_p.*t) + x1(2*p,1).*sin(omega_p.*t);
    v1_t = v1_t + omega_p*(x1(2*p,1).*cos(omega_p.*t) - x1(2*p-1,1).*sin(omega_p.*t));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Friction force, displacement, velocity in frequency domain: averaged
% harmonics cn^2=an^2+bn^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nphi_plot = 5; % Number of harmonics shown (not necessarily equal to Nphi)
% Nphi_plot should be smaller or equal to NPhi
rh=zeros(Nphi_plot,1); x2h=zeros(Nphi_plot,1); v2h=zeros(Nphi_plot,1);
x1h=zeros(Nphi_plot,1); v1h=zeros(Nphi_plot,1);

for p = 1:Nphi_plot
    omega_p = (2*p-1)*omega;

    % friction force
    rh(p) = sqrt(r(2*p-1,1)^2 + r(2*p,1)^2);

    % dof 2
    x2h(p) = sqrt(x2(2*p-1,1)^2 + x2(2*p,1)^2);
    v2h(p) = omega_p*sqrt(x2(2*p-1,1)^2 + x2(2*p,1)^2);
    
    % dof 1
    x1h(p) = sqrt(x1(2*p-1,1)^2 + x1(2*p,1)^2);
    v1h(p) = omega_p*sqrt(x1(2*p-1,1)^2 + x1(2*p,1)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modal projection of time histories and harmonics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[P,D] = eig(K,M);
disp = zeros(2,length(t));
disp(1,:) = x1_t;
disp(2,:) = x2_t;
disp_modal = P*disp;

velo = zeros(2,length(t));
velo(1,:) = v1_t;
velo(2,:) = v2_t;
velo_modal = P*velo;

% FFT of disp_modal
vec = 2*fft(disp_modal(1,1:end-1))/nt; %FFT to get the frequency domain
x1_modal = zeros(2*Nphi,1);
x1_modal(1:2:end-1,1) = real(vec(2:2:2*Nphi));
x1_modal(2:2:end,1) = -imag(vec(2:2:2*Nphi));

vec = 2*fft(disp_modal(2,1:end-1))/nt; %FFT to get the frequency domain
x2_modal = zeros(2*Nphi,1);
x2_modal(1:2:end-1,1) = real(vec(2:2:2*Nphi));
x2_modal(2:2:end,1) = -imag(vec(2:2:2*Nphi));

vec = 2*fft(velo_modal(1,1:end-1))/nt; %FFT to get the frequency domain
v1_modal = zeros(2*Nphi,1);
v1_modal(1:2:end-1,1) = real(vec(2:2:2*Nphi));
v1_modal(2:2:end,1) = -imag(vec(2:2:2*Nphi));

vec = 2*fft(velo_modal(2,1:end-1))/nt; %FFT to get the frequency domain
v2_modal = zeros(2*Nphi,1);
v2_modal(1:2:end-1,1) = real(vec(2:2:2*Nphi));
v2_modal(2:2:end,1) = -imag(vec(2:2:2*Nphi));

x1h_modal=zeros(Nphi_plot,1);
x2h_modal=zeros(Nphi_plot,1);
v1h_modal=zeros(Nphi_plot,1);
v2h_modal=zeros(Nphi_plot,1);
for p = 1:Nphi_plot
    omega_p = (2*p-1)*omega;
    
    x1h_modal(p) = sqrt(x1_modal(2*p-1,1)^2+x1_modal(2*p,1)^2);
    x2h_modal(p) = sqrt(x2_modal(2*p-1,1)^2+x2_modal(2*p,1)^2);
    v1h_modal(p) = sqrt(v1_modal(2*p-1,1)^2+v1_modal(2*p,1)^2);
    v2h_modal(p) = sqrt(v2_modal(2*p-1,1)^2+v2_modal(2*p,1)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% main title
sgtitle(['N_{\phi}=', num2str(Nphi), ', \mu = ', num2str(mu), ', N = ', num2str(N), ', \omega = ', num2str(omega)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Top row: time-domain histories over one period of forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,3,1)
plot(t,x1_t,'b',t,x2_t,'k')
ylabel('displacement (m)')
title('Physical displacements')
legend('x1','x2','Location','southeast')
xlabel('time (s)')
axis tight

subplot(5,3,2)
plot(t,v1_t,'b',t,v2_t,'k')
title('Physical velocities')
ylabel('velocity (m/s)')
legend('v1','v2','Location','southeast')
xlabel('time (s)')
axis tight

subplot(5,3,3)
plot(t,fric,'r')
title('Friction force')
xlabel('time (s)')
ylabel('r (N)')
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second row: physical response harmonics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,3,4);
hold on
bar((1:2:2*Nphi_plot)-0.15,x1h,'b','barwidth',0.15)
bar((1:2:2*Nphi_plot)+0.15,x2h,'k','barwidth',0.15)
%title('Physical displacements')
ylabel('participation')
xlabel('harmonic index')
xticks((1:2:2*Nphi_plot))
box on
axis tight

subplot(5,3,5)
hold on
bar((1:2:2*Nphi_plot)-0.15,v1h,'b','barwidth',0.15)
bar((1:2:2*Nphi_plot)+0.15,v2h,'k','barwidth',0.15)
%title('Physical velocities')
ylabel('pariticipation')
xlabel('harmonic index')
xticks((1:2:2*Nphi_plot))
box on
axis tight

subplot(5,3,6)
bar((1:2:2*Nphi_plot),rh,'r','barwidth',0.15);
%title('Friction force')
ylabel('participation')
xlabel('harmonic index')
xticks((1:2:2*Nphi_plot))
box on
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third row: time-domain histories in modal space 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,3,7)
plot(t,disp_modal(1,:),'b',t,disp_modal(2,:),'k')
title('Modal displacements')
ylabel('participation')
xlabel('time (s)')
legend('mode 1','mode 2','Location','southeast')
axis tight

subplot(5,3,8)
plot(t,velo_modal(1,:),'b',t,velo_modal(2,:),'k')
title('Modal velocities')
ylabel('participation')
xlabel('time (s)')
legend('mode 1','mode 2','Location','southeast')
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth row: modal response harmonics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,3,10)
hold on
bar((1:2:2*Nphi_plot)-0.15,x1h_modal,'b','barwidth',0.15);
bar((1:2:2*Nphi_plot)+0.15,x2h_modal,'k','barwidth',0.15);
%title('Modal displacements')
ylabel('participation')
xlabel('harmonic index')
xticks((1:2:2*Nphi_plot))
box on
axis tight

subplot(5,3,11)
hold on
bar((1:2:2*Nphi_plot)-0.15,v1h_modal,'b','barwidth',0.15);
bar((1:2:2*Nphi_plot)+0.15,v2h_modal,'k','barwidth',0.15);
%title('Modal velocities')
ylabel('participation')
xlabel('harmonic index')
xticks((1:2:2*Nphi_plot))
box on
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom row: vairous diagrams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,3,13)
plot(v2_t,fric,'k')
title("Coulomb's set")
ylabel('r(t)')
xlabel('v2(t)')
axis tight

subplot(5,3,14)
plot(x1_t,fric,'k')
title('Hysteresis 1')
ylabel('r(t)')
xlabel('x1(t)')
axis tight

subplot(5,3,15)
plot(x2_t,fric,'k')
title('Hysteresis 2')
ylabel('r(t)')
xlabel('x2(t)')
axis tight
