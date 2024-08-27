close all
clear
clc

rho_w = 1026; % density of liquid (kg/m^3) [water]
P_atm = 101.325e3; % atmospheric pressure
gamma = 1.299; % heat ratio

R0 = [8e-5 1e-4 3e-3 0.0446]; % bubble radius (m)
f0 = [30e3 70e3 72.8e3 75e3]; % frequency range 30k to 70k (Hz)

% Minnaert frequency (resonance frequency)
f_res = 1./(2*pi*R0)*sqrt(3*gamma*P_atm/rho_w)
R0_res = 1./(2*pi*f0)*sqrt(3*gamma*P_atm/rho_w)


% function f_res = minnaert_freq1(R0)
%     rho_w = 1026; % density of liquid (kg/m^3) [water]
%     P_atm = 101.325e3; % atmospheric pressure
%     gamma = 1.299; % heat ratio
%     f_res = 1./(2*pi*R0)*sqrt(3*gamma*P_atm/rho_w);
% end
% 
% function R0 = minnaert_radius(f_0)
%     rho_w = 1026; % density of liquid (kg/m^3) [water]
%     P_atm = 101.325e3; % atmospheric pressure
%     gamma = 1.299; % heat ratio
%     R0_res = 1./(2*pi*f0)*sqrt(3*gamma*P_atm/rho_w);
% end