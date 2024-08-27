% Date: 02.04.2024
% Name: sonar equation implementation
% Wave propagation simulation from sonar to the bubble with the TS calculation
clc; clear; close all

% SL = 220; 
P_ref = 10^5; % Pa, power at the ref. distance of 1m
SL = 10*log10(P_ref*0.08/(0.67*1e-18)); % dB re 1mPa at 1m, source level
DI = 10*log10(32); % dB, directivity index is the gain in SNR 
AG = DI; % dB, array gain
R = 10; % m, range
alpha_att = 0.11; % dB/km, absorption coefficient

%% NL noise level
NL = 50; % dB, noise level
T = 0.1; % s, pulse length
BW = 1/T; % bandwidth of the receiver
NL = NL + 10*log10(BW);

%% TS target strength

f = 75000; % Hz=1/s , frequency
c = 1500; % m/s, speed of sound in water
r0 = 4e-3; % bubble radius
% r0 = 3e-4;

sigma_bs = bubble_response_model(f,r0,1);
% TS_total = - 2*TL + TS 
TS_t =  TS_total(R, alpha_att, sigma_bs) ; % dB, total target strength
TS2 = 10*log10(r0^2/4); % dB, target strength


SNR = SL + TS_t - (NL - AG) % signal to noise ratio, active sonar equation 
%% 0. General the transmission loss
function TL = TLsonar_gen(R, alpha_att)
    TLg = 20*log10(R);% spherical TL from geometrical spreading, one way
    TLl = alpha_att * (R); % dB, attenuation ~ TL from dissipation
    TL = TLg + TLl; % dB, transmission loss 
end
%% 1. Calculate the transmission loss from SONAR to 1m in front of target
function TL = TLsonar_to1m_tar(R, alpha_att)
    TLg = 20*log10(R-1);% spherical TL from geometrical spreading, one way
    TLl = alpha_att * (R-1); % dB, attenuation ~ TL from dissipation
    TL = TLg + TLl; % dB, transmission loss 
end
%% 2. Calculate TL from 1m in front of target to target
function TL = TLtar_to1m_tar(alpha_att)
    TLg = 20*log10(1);% spherical TL from geometrical spreading, one way
    TLl = alpha_att * (1); % dB, attenuation ~ TL from dissipation
    TL = TLg + TLl; % dB, transmission loss 
end
%% 3. Calculate TS from filtering with frequency response of target
function TS = TS_tar(sigma_bs)
    TS = min(10*log10(sigma_bs)); 
end
%% 4. Add transmission loss from target to 1 m in front of target
% - Combine these for total target strength
function TS_t = TS_total(R, alpha_att, sigma_bs) 
 % TS_total = - 2*TL + TS
    TL_t0 = -(TLsonar_to1m_tar(R, alpha_att) + TLtar_to1m_tar(alpha_att)) * 2
    TL_t = -TLsonar_gen(R, alpha_att) * 2
    TS = TS_tar(sigma_bs)
    TS_t = TL_t0 + TS; 
end
