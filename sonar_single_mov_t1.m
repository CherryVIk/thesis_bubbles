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
%% Signal Parameters
% Sample frequency
fs = 192000;
% Signal bandwidth
fB = 40000;
% Center frequency
fC = 50000;
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
    tx = filter(Hd, tx);
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
fig=figure;
subplot(211);
grid on;
plot(t, tx(:,1));
xlabel('Time, s');ylabel('Amplitude');
title("Transmitted, time domain signal")
best_plot_ever(fig)

subplot(212);
logTx = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
grid on;
plot(-f(1:NBins), logTx)
xlabel('Frequency, Hz');ylabel('Magnitude (dB)');
title("Transmitted, log spectrum ")
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov_transmitted_50kHz_time+freq","png");

%% Environment settings

% zcoord = -50 + move_ii;
zcoord = 0;
% posTar = [0 10 zcoord; 20 10 0]; 
posTar = [0 0 20];
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
noise_level_dB = -65;
noise_level_linear = 10^(noise_level_dB/10);
noise_add = randn(nRxSeqLength, NRx) * noise_level_linear; 
% rx = rx + noise_add;

radius_b = 8e-5;%3e-3;%8e-5;%3e-4;% Oscillations, bubble radius (m)
str_radius = num2str(radius_b);
str_radius = strrep(str_radius, '.', ',');
sigma_bs = bubble_response_model(f,radius_b,1);

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
            rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) + mixed_resp(1:nSig, iTx);
  
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
%% Apply damping: 1/r (linear), 1/r^2 (cylindrical), 1/r^3 (spherical)
propagationTimesPerSample = (0:nRxSeqLength)./ fs;
propagationDistancePerSample = propagationTimesPerSample * cWater;
damping = 1./propagationDistancePerSample.^geoSpreadLoss;

%% Test correlation
% Calculate crosscorrelation between rx & tx signals
% For mid rx-element
[corr, lag] = xcorr(rx(:, round(NRx/2)), tx(:, round(NTx/2)));
tLag = lag;
% Normalize correlation to max
corr = abs(corr(nRxSeqLength:end) ./ max(abs(corr(nRxSeqLength:end))));
% Only consider positive correlation lags (total length = 2xnRxSeqLength)
tLagInMeters = tLag(nRxSeqLength:end)./fs.*cWater;
% Find max correlation peaks
[pk, loc] = findpeaks(corr, 'MinPeakHeight', 0.8, 'MinPeakDistance',6);
% Find shortest propagation time for mid tx & rx element
locShortest = min(squeeze(tPropagationTime(round(NTx/2), :, round(NRx/2))));
% Calculate time vector
tSim = linspace(0, nRxSeqLength/fs, nRxSeqLength);

% Plot test correlation --------------------------------------------------------------------

fig=figure;
subplot(211);
plot(tSim, rx(:, 1));
grid on;
title('Received signal');
xlabel("Time, s"); ylabel("Amplitude")
best_plot_ever(fig)
subplot(212);
semilogx(tLagInMeters, corr);
hold on;
semilogx(tLagInMeters(loc), pk, 'rx');
grid on;
xlabel("Distance, m"); 
title('Crosscorrelation: Transmit- $\&$ receive signal');
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov+noise_correlated_50kHz_r"+str_radius,"png");

%% Test Bubble Response Plot --------------------------------------------------------------------
f_sigma_bs = abs(Tx(:, iTx)).*sigma_bs(1:NBins);
res_temp = f_sigma_bs./abs(Tx(:, iTx));

fig=figure;
logFs = 10*log10(abs(f_sigma_bs(:,1))./max(abs(f_sigma_bs(:,1))));
hold on
plot(-f(1:NBins),logFs)
logRes_s = 10*log10(res_temp(1:NBins,1));
plot(-f(1:NBins), logRes_s,'-.','LineWidth',1.5);
xlim([0 f(end)])
ylim([-100 0])
title("Transmitted signal with bubble response, log spectrum");
subtitle("r = " + str_radius + " m")
xlabel("Frequency, Hz");ylabel("Magnitude, dB")
legend('transmitted signal with bubble response' ,'bubble response','Location', 'Best')

grid on;
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov+transmitWithBubbleResponse-set_50kHz_r"+str_radius,"png");
%% Plot single bubble response (time and freq domain)
Y = fft(rx, NFFT);
Y = Y(1:NBins, :);

fig=figure;
logRx_withB = 20*log10(abs(Y(:,1))./max(abs(Y(:,1))));
plot(-f(1:NBins), logRx_withB(:, 1));
grid on;
hold on
xlabel('Frequency (Hz)');ylabel("Magnitude, dB")
title('Received signal, spectrum');
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov_received_50kHz_spectrum_noise_r"+str_radius,"png");

%% Reconstruction of the signal of before the correlation
Y = Y.*hamming(length(Y));
H_hat = abs(Y(:,1)) ./ abs(Tx(:,1));

fig=figure;
f_half =  -f(1:NBins);

hold on
logH_hat = 20*log10(abs(H_hat(:,1))./max(abs(H_hat(:,1))));
[pk, loc] = findpeaks(logH_hat,'MinPeakHeight', -1e-4);
plot(-f(1:NBins), logH_hat(:, 1));
plot(-f(1:NBins),logRes_s)
hold on;
semilogx(f_half(loc), pk,'rx');
xline(f_half(loc),'r');
ylim([-100 0])
legend("reconstruction","bubble")
xlabel('Frequency (Hz)');ylabel("Magnitude, dB")
title('Reconstructed spectrum');
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov+noise_reconstructed_with-peak_r"+num2str(radius_b),"png");

%% Beamforming: Calculate array manifold vector (AMV)
NFFT = 2^nextpow2(2 * nRxSeqLength); 
% Frequency support vector: freq = bin*fs/FFT_size
NBins = NFFT / 2 + 1; % FFT-bins of pos. frequencies
bins = 0:NBins-1; % Freq. bin support vector
f = bins * fBinRes;
% Calculate delays wrt. array center for all beamforming angles
tTau = (sind(angles) .* posRx(:,1)) / cWater;
% Create array manifold vector [iNumEl,iFftHalfPlusOne,iNumBeams]
% Define the dimensions of the input arrays

% Compute the array manifold vector
% Repeat the frequency vector for all [beams, freqs, NRx]
Frame=repmat(f,[NBeams, 1,NRx]);
% Reshape to [beams, NRx, freqs]
Frame = permute(Frame, [1 3 2]);
% Multiply tau with freq. supports
tauTimesf = bsxfun(@times,tTau', Frame);
% Create exponential terms: e^(j2pi*tau*f)
AMV = exp(1j * 2 * pi *tauTimesf);
AMV = permute(AMV, [2 3 1]);

%% Matched filtering: Correlate receive & transmit signals
% rx signal to freq. domain
Rx = fft(rx, NFFT);
Rx = Rx(1:NBins, :);

% tx signal to freq. domain
Tx = fft(tx, NFFT);
Tx = Tx(1:NBins, :);
% Matched filtering: Find tx-sequence in each rx-sequence
Mf = conj(Tx) .* Rx;

% Delay-and-sum beamforming: Combine all correlated rx-sequences by
% multiplication with the array manifold vector to focus the energy in the
% beam angle direction for each beam over -90° to 90°
MfBf = squeeze(sum(Mf.' .* AMV, 1)).';

% Transform beamformed & correlated output to time domain
mfbf = ifft(MfBf,NFFT, 2, 'symmetric');
mfbf = abs(mfbf(:, 1:nRxSeqLength));

% Normalize correlated output
mfbf = mfbf ./ max(max(mfbf));
NResized = 1000;
%% Downscaling: Summarize correlated cells by their max values to 
% reduce overall output size
% Downscaling factor: Summarize N distance cells by their max value
downScale = nRxSeqLength / NResized;
mfbfsmall = zeros(NBeams, NResized);
for iResize = NResized:-1:1
    startIndex = round((iResize-1)*downScale + 1);
    endIndex = round(iResize*downScale);
    mfbfsmall(:,iResize) = ...
        max( ...
        mfbf(:,startIndex:endIndex ...
        ),[],2);
end

%% Reconstruction of the correlated signal
Y = Mf;
Y = Y.*hamming(length(Y));
H_hat = abs(Y(:,1)) ./ abs(Tx(:,1));
% Normalize to max
logH_hat = 20*log10(abs(H_hat(:,1))./max(abs(H_hat(:,1))));
% H_hat = H_hat.*hamming(length(H_hat));

f = linspace(-fs/2, fs/2, NFFT);
f_half =  -f(1:NBins);
%% Identify the bubble resonance peak
% Find max peaks
[pk, loc] = findpeaks(logH_hat,'MinPeakHeight', -1e-4);

tSim = linspace(0, nRxSeqLength/fs, nRxSeqLength);
% Plot --------------------------------------------------------------------
fig=figure;
subplot(211);
plot(tSim, rx(:, 1));
grid on;
title('Received signal');
xlabel("Time, s");ylabel("Amplitude");
best_plot_ever(fig)
subplot(212);
plot(f_half, logH_hat);
% ylim([0 3e-4])
hold on;
semilogx(f_half(loc), pk,'rx');
xline(f_half(loc),'r')
grid on;
xlabel("Frequency, Hz"); ylabel("Magnitude, dB"); 
title('Reconstructed correlated signal with a peak');
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov_received_50kHz_with-peak","png");
%% Find bubble properties 
% the resonance frequency
f_0 = f_half(loc);
[pk, loc2] = findpeaks(logRes_s,'MinPeakHeight', -63)
f_0_input = f_half(loc2);

f0_minn = minnaert_freq(radius_b);
disp("input: f0 = " + num2str(f0_minn) + newline + " and f0_input = " + num2str(f_0_input)  + newline + "output f0 = " + num2str(f_0));

% the bubble radius
R0 = minnaert_radius(f_0);
disp("input: R0 = " + str_radius + newline + "output R0 = " + num2str(R0));
%% Plot the reconstructed correlated signal
% f = linspace(-fs/2, fs/2, NFFT);
fig=figure;
subplot(311)
% logRx_noB = X(:,1);
logTx_noB = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
plot(-f(1:NBins), logTx_noB(:, 1));
grid on;
title('Transmitted signal, freq');
best_plot_ever(fig)

subplot(312)
% logRx_withB = Y(:,1);
logRx_withB = 20*log10(abs(Y(:,1))./max(abs(Y(:,1))));
plot(-f(1:NBins), logRx_withB(:, 1));
grid on;
hold on
title('Received correlated signal, freq');
best_plot_ever(fig)

subplot(313)
hold on
logH_hat = 20*log10(abs(H_hat(:,1))./max(abs(H_hat(:,1))));
plot(-f(1:NBins), logH_hat);
plot(-f(1:NBins), logRes_s, '-.', 'LineWidth', 1.5);
ylim([-100 0])
legend("reconstruction", "bubble")
title('Reconstructed correlated signal, freq');
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov_reconstr_corr_with-peak_set_50kHz_r"+str_radius,"png");
%% 
fig=figure;
f_half =  -f(1:NBins);

hold on
logH_hat = 20*log10(abs(H_hat(:,1))./max(abs(H_hat(:,1))));
[pk, loc] = findpeaks(logH_hat,'MinPeakHeight', -3e-1);
plot(-f(1:NBins), logH_hat(:, 1));
plot(-f(1:end/2+1),logRes_s, '-.', 'LineWidth', 1.5)
hold on;
semilogx(f_half(loc(4)), pk(4),'rx');
xline(f_half(loc(4)),'r');
ylim([-100 0])
legend("reconstruction","bubble",'Location','best')
title('Reconstructed correlated signal, freq');
xlabel("Frequency, Hz"); ylabel("Magnitude, dB"); 
best_plot_ever(fig)
% saveas(gca, "thesis_pics/single_mov_reconstr_corr_with-peak_50kHz_r"+str_radius,"png");

%% Plot results as PPI (plan position indicator) plot
% Log of correlated output
sonar_fig=figure;
ppi = 20*log10(abs(mfbfsmall)+eps);
ppi(ppi<-40) = -40;
% Maximum distance (single path) = Half of max sequence travel distance
dMax = nRxSeqLength / fs * cWater / 2; 
theta = deg2rad(angles + 90);
r = 0:dMax/size(ppi,2):dMax-1/size(ppi,2);
[THETA,RR] = meshgrid(theta,r);
[A,B] = pol2cart(THETA,RR);
surf(A,B,ppi.');
colormap('jet');
view(0,90);
xlabel('x [m]');
ylabel('y [m]');
daspect([1 1 1]);
axis tight
shading interp;
colorbar;
ax = gca;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'units', 'pixels', 'position', [100 40 1500 900]);
hold on;
hold off;
best_plot_ever(sonar_fig)
set(gca,'FontSize',12)
% saveas(gca, "thesis_pics/single_mov_sonarFig_50kHz_r"+num2str(radius_b),"png");

%% Functions
function f_res = minnaert_freq(R0)
    rho_w = 1026; % density of liquid (kg/m^3) [water]
    P_atm = 101.325e3; % atmospheric pressure
    gamma = 1.299; % heat ratio
    f_res = 1./(2*pi*R0)*sqrt(3*gamma*P_atm/rho_w);
end

function R0_res = minnaert_radius(f_0)
    rho_w = 1026; % density of liquid (kg/m^3) [water]
    P_atm = 101.325e3; % atmospheric pressure
    gamma = 1.299; % heat ratio
    R0_res = 1./(2*pi*f_0)*sqrt(3*gamma*P_atm/rho_w);
end