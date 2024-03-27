% plotting script. Plots over the last period assuming steady state is reached.
% Is called by IntTime_Nonlin.m

nTT = floor(nT/div);
figure
sgtitle('Nonlinear simulation')

subplot(2,2,1)
plot(time(1:nTT+1),sol(1,end-nTT:end),'b',time(1:nTT+1),sol(3,end-nTT:end),'k')
xlabel('time')
legend('x1', 'x2','Location','southeast')
ylabel('displacement')
axis tight

subplot(2,2,2)
plot(time(1:nTT+1),sol(2,end-nTT:end),'b',time(1:nTT+1),sol(4,end-nTT:end),'k')
xlabel('time')
legend("x1'","x2'",'Location','southeast')
ylabel('velocity')
axis tight

subplot(2,2,3)
plot(time(1:nTT+1),sol(5,end-nTT:end),'r')
xlabel('time')
ylabel('friction force r')
axis tight

subplot(2,2,4)
% note that because of the storage procedure, Coulomb's set 
% might have a strange look in the plot. In "IntTime_PredCorr.m, set
% "div=1" to recover the correct look.
plot(sol(4,end-nTT:end),sol(5,end-nTT:end),'k')
xlabel("velocity x2'")
ylabel('friction force r')
axis tight