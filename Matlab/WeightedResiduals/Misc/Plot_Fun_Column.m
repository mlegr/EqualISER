% Solution plot over one period
% Organises the plots wolumnwise

% function plots [works]
t = linspace(0,2*pi/omega,2000);

fric = zeros(size(t));

x2_t = zeros(size(t));
v2_t = zeros(size(t));

x1_t = zeros(size(t));
v1_t = zeros(size(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friction force and displacement
% harmonics (cos/sin)
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

figure
subplot(4,1,1)
plot(t,x2_t,'k',t,x1_t,'b')
ylabel('displacement')
%legend('x2','x1','Location','southeast')
xlabel('time')
axis tight

subplot(4,1,2)
%plot(t,v2,'k',t,fric,'r')
plot(t,v2_t,'k',t,v1_t,'b')
%ylabel('norm. mag.')
ylabel('velocity')
%legend('v2','v1','Location','southeast')
xlabel('time')
axis tight

subplot(4,1,3)
%plot(x2,v2,'k')
plot(t,fric,'r')
%title('phase diagram')
xlabel('t')
ylabel('r')
axis tight
%ylim([min(v2)-0.2 max(v2)+0.2])
%xlim([min(x2)-0.2 max(x2)+0.2])

subplot(4,1,4)
plot(v2_t,fric,'k')
ylabel('r(t)')
xlabel('v2(t)')
axis tight

 set(gcf, 'Position',  [1500, 100, 250, 1000])

%%%%%%%%%%%%%%%%%%%%%%%%
%figure(2)
%subplot(2,1,1)
%plot(t,x1,'r',t,x2,'k')
%title('displacements')
%legend('dof1','dof2')
%xlabel('time')
%ylabel('displacement')
%axis tight
%ylim([-0.25 0.25])

%subplot(2,1,2)
%plot(t,v1,'r',t,v2,'k')
%title('velocities')
%legend('dof1','dof2')
%xlabel('time')
%ylabel('velocity')
%axis tight
%ylim([-1.5 1.5])
%filename = ['omega' num2str(omega) 'dofNphi' num2str(Nphi) '.pdf'];
%saveas(2,filename)