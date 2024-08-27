clear all;
close all;
%% SONAR Parameters
% Speed of sound underwater
cWater = 1500.0;
% Number of transmitting projectors
NTx = 1; % SIMO: Single-input
centerTx = [0 0 0];
% Number of receiving hydrophones
NRx = 32; % SIMO: Multi-output
centerRx = [0 0 0];
% Beamforming angles
angles = -90:2:90;
% Number of beams
NBeams = length(angles);

filename = 'MyAnimation.gif';
% for move_ii = 100
%% Signal Parameters
% Sample frequency
fs = 192000;
% Signal bandwidth
fB = 20000;
% Center frequency
fC = 75000;
% Min & max frequencies
fMin = fC - fB/2;
fMax = fC + fB/2;
% Ping duration
tPing = 1.0; % in seconds
nPing = tPing * fs; % in samples
% Signal duration
tSig = 0.05; % in seconds
nSig = tSig * fs; % in samples
t = linspace(0, tSig, nSig);
% Signal types
eSignalTypes.CW = 'CW';
eSignalTypes.blNoise = 'blNoise';
eSignalTypes.HFM = 'HFM';
eSignalType = eSignalTypes.blNoise;
bRandomSeed = 42;%42
rng(bRandomSeed)

%% Adaption of SONAR-system wrt. signal design
% Set array element distance to half of wavelength of max freq.
dRx = cWater / (fMax) / 2; % Lambda half distance of array elements
dTx = dRx;

%% FFT-parameters
nNextPow2 = nextpow2(nSig*2); % find nearest x so 2^x = nSig
NFFT = 2^nNextPow2; % FFT-length as multiple of 2
% NFFT = nSig;
NBins = NFFT / 2 + 1; % FFT-bins of pos. frequencies
bins = 0:NBins-1; % Freq. bin support vector
fBinRes= fs / NFFT;
nfMin = floor(fMin/fBinRes);
nfMax = ceil(fMax/fBinRes);
f = linspace(-fs/2, fs/2, NFFT);%NFFT
% f = linspace(0, fs/2, NBins);%NFFT

%% Generate transmit sequence
% Bandpass filter design
Fstop1 = fMin-100;       % First Stopband Frequency
Fpass1 = fMin;       % First Passband Frequency
Fpass2 = fMax;       % Second Passband Frequency
Fstop2 = fMax+100;       % Second Stopband Frequency
Astop1 = 100;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 100;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, fs);
Hd = design(h, 'cheby2');



if strcmp(eSignalType, eSignalTypes.blNoise)
    % Generate Gaussian white noise
    tx = randn(nSig, NTx);
    % an amplitude for that noise is 10% of the noise-free signal at every element.
     % noise_add = 
   
    tx = filter(Hd, tx);
    % tx = tx + noise_add;
    %tx = filtfilt(Hd.sosMatrix, Hd.ScaleValues, tx);
    % Transform time to freq. domain signal
    Tx = fft(tx, NFFT);%NFFT
    % Only save positive freq.
    Tx = Tx(1:NBins, :);
% The following has a commented out ideal (but impractical) bandpass filter
    % Bandpass
    %Tx(1:nfMin, :) = 0;
    %Tx(nfMax:end, :) = 0;
    % Freq. to time domain
    %tx = ifft(Tx, NFFT, 'symmetric');
    %tx = tx(1:nSig);
% End of ideal bandpass example
end

%% Plot transmit sequence
figure;
subplot(211);
plot(t, tx(:,1));
title("Time domain signal");
grid on;
subplot(212);
logTx = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
plot(f(end/2:end), logTx);
title("Log. frequency spectrum ");
grid on;

%% Environment settings

% zcoord = -50 + move_ii;
zcoord = 0;
% posTar = [0 10 zcoord; 20 10 0]; 
posTar = [0 40 0];
NTargets = size(posTar, 1);
bDirectSound = 0;

%% SONAR-sytem element positions (line arrays)
% Uniform line array element equidistant element spacing 
% around array center (spaced on x-axis)
rxCenterElement = (NRx + 1) / 2;
txCenterElement = (NTx + 1) / 2;
turnRx = [cosd(0) sind(0) 0];
turnTx = [cosd(0) sind(0) 0];
posRx = (((1:NRx) - rxCenterElement) .* dRx)' .*turnRx+ centerRx;
posTx = (((1:NTx) - txCenterElement) .* dTx)' .*turnTx+ centerTx;

%% Calculate propagation delays
tPropagationDist = zeros(NTx, NTargets+bDirectSound, NRx);
% Propagation distance: Tx -> target -> Rx
for iTx = 1:NTx
    for iRx = 1:NRx
        for iTar = 1:NTargets
        dTxToTar = norm(posTar(iTar, :) - posTx(iTx, :));
        dTarToRx = norm(posRx(iRx, :) - posTar(iTar, :));
        tPropagationDist(iTx, iTar, iRx) = dTxToTar + dTarToRx;
        end
        if (bDirectSound == 1)
        tPropagationDist(iTx, iTar + bDirectSound, iRx) = ...
        norm(posRx(iRx, :) - posTx(iTx, :));
        end
    end
end

% Propagation distance -> Propagation time
tPropagationTime = tPropagationDist ./ cWater; % in seconds
tPropagationTime = tPropagationTime .* fs; % in samples

%% Calculate received signals
% Geometric spreading loss damping
geoSpreadLoss = 0;
% geoSpreadLoss = 1/( norm(posTar(iTar, :) - posTx(iTx, :))^2 ); % 1/r^x power loss
% Max rx. sequence length (signal duration + max propagation time)
nRxSeqLength = nSig + ceil(max(tPropagationTime(:)));
rx = zeros(nRxSeqLength, NRx);

%% Background noise
noise_level_dB = -60;
noise_level_linear = 10^(noise_level_dB/10);
noise_add = randn(nRxSeqLength, NRx) * noise_level_linear; 
rx = rx + noise_add;

radius_b = 2e-4;% 2e-4, 2e-3 Oscillations, bubble radius (m)
sigma_bs = bubble_response_model(f,radius_b, 1);
% 
% Plot --------------------------------------------------------------------
f_sigma_bs = abs(Tx(:, iTx)).*sigma_bs(1:NBins);
res_temp = f_sigma_bs./abs(Tx(:, iTx));

figure;
subplot(211);
logBs = 10*log10(sigma_bs(1:NBins,1));
plot(f(end/2:end), logBs);
ylim([-200 0])
hold on
logRes_s = 10*log10(res_temp(1:NBins,1));
plot(f(end/2:end), logRes_s);
title("Log. frequency spectrum, bubble");
grid on;
subplot(212);
logFs = 10*log10(abs(f_sigma_bs(:,1))./max(abs(f_sigma_bs(:,1))));
plot(f(end/2:end),logFs)
ylim([-200 0])
title("Log. frequency spectrum, with bubble ");
grid on;
%%
% sigma_bs=ones(NBins,NTx);
for iTx = 1:NTx
    for iRx = 1:NRx
        for iTar = 1:NTargets + bDirectSound
            iStart = floor(tPropagationTime(iTx, iTar, iRx))+1;
            iEnd = floor(iStart + length(tx(:, iTx))) - 1;
            % Add reflected tx-signal at delay time to rx_-
%             rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) +  tx(:, iTx);
            f_sigma_bs = abs(Tx(:, iTx)).*sigma_bs(1:NBins);
            theta = angle(Tx(:, iTx));
            mixed_resp = f_sigma_bs.*exp(i*theta);
            mixed_resp = ifft(mixed_resp, NFFT,'symmetric');
            mixed_resp = mixed_resp(1:nSig);
%             idea: to add phase shift (imag part of the Tx to the mixed
%             freq resp)
%           rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) + tx(1:nSig, iTx);
            rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) + mixed_resp(1:nSig, iTx);

            

%% Insert custom freq. response here!   
% As of now, the transmission signal is simply added to the receive signal
% but as the reflection at the object happens, in the time domain, the transmit sequence is
% convolved with the impulse response of the target. In the
% frequency domain, the spectrum of the transmission signal is multiplied
% with the frequency response of the target.

%  filter your Tx signal with this frequency response for a given bubble radius

%% Insert custom freq. response here!   
        end
    end
end

%% Reconstruction of the signal
Y = fft(rx, NFFT);
Y = Y(1:NBins, :);

% [~,ff_s] = min(abs(f-Fstop1));
% [~,ff_e] = min(abs(f-Fstop2));
% ff_s = ff_s - NBins;
% ff_e = ff_e - NBins;
% ff_s = 10000;
% ff_e = 15000;
% H_hat = R / T
H_hat = abs(Y(:,1)) ./ abs(Tx(:,1));
% H_hat = abs(H_hat);


figure;
subplot(311)
% logRx_noB = X(:,1);
logTx_noB = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
plot(f(end/2:end), logTx_noB(:, 1));
grid on;
title('Target signal, freq');
subplot(312)
% logRx_withB = Y(:,1);
logRx_withB = 20*log10(abs(Y(:,1))./max(abs(Y(:,1))));
plot(f(end/2:end), logRx_withB(:, 1));
grid on;
hold on
title('Target-bubble signal, freq');
subplot(313)
% logH_hat = H_hat(:,1);
logH_hat = 20*log10(abs(H_hat(:,1)./max(abs(H_hat(:,1)))));
plot(f(end/2:end), logH_hat(:, 1));
hold on
logRes_s = 10*log10(res_temp(1:NBins,1));
plot(f(end/2:end), logRes_s);
title('Extracted bubble signal, freq');

%%
figure;
subplot(311)
plot(f(end/2:end), abs(Tx(:, 1)));
grid on;
title('Target signal, no log, freq');
subplot(312)
logRx_withB = 20*log10(abs(Y(:,1))./max(abs(Y(:,1))));
plot(f(end/2:end), Y(:, 1));
grid on;
hold on
title('Target-bubble signal, no log, freq');
subplot(313)
logH_hat = H_hat(:,1);
logH_hat = 20*log10(abs(H_hat(:,1))./max(abs(H_hat(:,1))));
plot(f(end/2:end), H_hat(:, 1));
hold on
plot(f(end/2:end), res_temp);
title('Extracted bubble signal, no log, freq');
%% Wiener filter
x = noise_add(:,1); y = rx(:,1); 
y = y(:); % reference signal
x = x(:); % signal with additive Gaussian noise
N = 200; % filter order
[xest,b,MSE] = wienerFilt(x,y,N);
%% plot results

tSim = linspace(0, nRxSeqLength/fs, nRxSeqLength);

figure
subplot(311)
plot(tSim(N+1:end),x(N+1:end),'r')
hold on
plot(tSim(N+1:end),y(N+1:end),'k')
ylim([-1e-5,1e-5]);
title('Wiener filtering example')
legend('noisy signal','reference')
subplot(312)
plot(tSim(N+1:end),xest,'k')
ylim([-1e-5,1e-5]);
legend('estimated signal')
subplot(313)
plot(tSim(N+1:end),(x(N+1:end) - xest),'k')
legend('residue signal')
xlabel('time (s)')
ylim([-1e-5,1e-5]);

%% Functions
% https://de.mathworks.com/matlabcentral/fileexchange/71440-signal-separation-with-wiener-filtering
   
function [xest,B,MSE] = wienerFilt(x,y,N)
    X = 1/N .* fft(x(1:N));
    Y = 1/N .* fft(y(1:N));
    X = X(:);
    Y = Y(:);
    Rxx = N .* real(ifft(X .* conj(X))); % Autocorrelation function
    Rxy = N .* real(ifft(X .* conj(Y))); % Crosscorrelation function
    Rxx = toeplitz(Rxx);
    Rxy = Rxy';
    B = Rxy / Rxx; B = B(:); % Wiener-Hopf eq. B = inv(Rxx) Rxy
    xest = fftfilt(B,x);
    xest = xest(N+1:end); % cut first N samples due to distorsion during filtering operation
    MSE = mean(y(N+1:end) - xest) .^2; % mean squared error
end