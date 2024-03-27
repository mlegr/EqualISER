% Description: HBM forced-response for 2-dof oscillator with Coulomb's friction on the
% second dof. Creates a Frequency Response plot of the system.
% Calls Main_Freq.m
% Author: Mathias Legrand
% Date: March 22, 2024

clear
close all

% Delta of the forcing frequency omega
delta_w = 0.01;

% Omega range for the frequency plot
omega_ini = 0.1;
omega_final = 4;

% Loads the mechanical parameters of the system
run("../parameters")

kk = 0; %counter
for omega = omega_ini:delta_w:omega_final
    
    kk = kk + 1

    % Calls (adapated) main file
    Main_Freq;

    % Harmonic averages
    for k = 1:Nphi
        x2h(k,kk) = sqrt(x2(2*k-1,1)^2 + x2(2*k,1)^2);
        x1h(k,kk) = sqrt(x1(2*k-1,1)^2 + x1(2*k,1)^2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Total energy = potential + kinetic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % kinetic energy
    Ec1 = 0; % mass 1
    Ec2 = 0; % mass 2
    for k = 1:Nphi
        Ec1 = Ec1 + (2*k-1)^2*(x1(2*k-1)^2+x1(2*k)^2);
        Ec2 = Ec2 + (2*k-1)^2*(x2(2*k-1)^2+x2(2*k)^2);
    end
    Ec1 = m1/4*Ec1;
    Ec2 = m2/4*Ec2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential energy
    Ep1 = 0; % spring 1
    Ep2 = 0; % spring 2
    for k = 1:Nphi
        Ep1 = Ep1 + x1(2*k-1)^2 + x1(2*k)^2;
        Ep2 = Ep2 + (x1(2*k-1)-x2(2*k-1))^2 + (x1(2*k)-x2(2*k))^2;
    end
    Ep1 = k1/4*Ep1;
    Ep2 = k2/4*Ep2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Total energy
    energy(kk,1) = Ec1 + Ec2 + Ep1 + Ep2;
end

figure
hold on
plot((omega_ini:delta_w:omega_final), energy);
set(gca, 'YScale', 'log');
title('Frequency-Energy plot')
ylabel('Energy [J]')
xlabel('Frequency [rad/s]')
box on