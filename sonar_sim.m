%%% NAME: Single bubble sonar simulation
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
fB = 20000;
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
bRandomSeed = 42;

%% Adaption of SONAR-system wrt. signal design
% Set array element distance to half of wavelength of max freq.
dRx = cWater / (fMax) / 2; % Lambda half distance of array elements
dTx = dRx;

%% FFT-parameters
nNextPow2 = nextpow2(nSig*2); % find nearest x so 2^x = nSig
NFFT = 2^nNextPow2; % FFT-length as multiple of 2
NBins = NFFT / 2 + 1; % FFT-bins of pos. frequencies
bins = 0:NBins-1; % Freq. bin support vector
fBinRes= fs / NFFT;
nfMin = floor(fMin/fBinRes);
nfMax = ceil(fMax/fBinRes);
f = linspace(-fs/2, fs/2, NFFT);%NFFT

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
%%
R = round(NFFT/nSig);
if strcmp(eSignalType, eSignalTypes.blNoise)
    tx = randn(nSig, NTx);
    tx = filter(Hd, tx);
    %tx = filtfilt(Hd.sosMatrix, Hd.ScaleValues, tx);
    % Transform time to freq. domain signal
    Tx = fft(tx, NFFT);
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

% Plot transmit sequence
figure;
subplot(211);
%t = linspace(0,tSig,NFFT);
plot(t, tx(:,1));
title("Time domain signal");
grid on;
subplot(212);
logTx = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
plot(f(end/2:end), logTx);
title("Log. frequency spectrum ");
grid on;

%% Environment settings
posTar = [0 40 0 ; -60 20 -10 ;10 20 0]; % [x y z]
NTargets = size(posTar, 1);
bDirectSound = 0;

%% SONAR-sytem element positions (line arrays)
% Uniform line array element equidistant element spacing 
% around array center (spaced on x-axis)
rxCenterElement = (NRx + 1) / 2;
txCenterElement = (NTx + 1) / 2;
turnRx = [cosd(0) sind(0) sind(0)];
turnTx = [cosd(0) sind(0) sind(0)];
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
geoSpreadLoss = 0; % 1/r^x power loss
% Max rx. sequence length (signal duration + max propagation time)
nRxSeqLength = nSig + ceil(max(tPropagationTime(:)));
rx = zeros(nRxSeqLength, NRx);

radius_b = 0.5e-3;% Oscillations, bubble radius (m)
sigma_bs = bubble_response(f,radius_b);
%%
f_sigma_bs = abs(Tx(:, iTx)).*sigma_bs(1:NBins);

figure;
subplot(211);
logBs = 10*log10(sigma_bs(1:NBins,1));
plot(f(end/2:end), logBs);
ylim([-100 0])
title("Log. frequency spectrum, bubble");
grid on;
subplot(212);
logFs = 20*log10(abs(f_sigma_bs(:,1))./max(abs(f_sigma_bs(:,1))));
plot(f(end/2:end),logFs)
ylim([-100 0])
title("frequency spectrum, with bubble ");
grid on;
%%
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
            rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) ;% + mixed_resp(1:nSig, iTx);

%% Insert custom freq. response here!   
% As of now, the transmission signal is simply added to the receive signal
% but as the reflection at the object happens, in the time domain, the transmit sequence is
% convolved with the impulse response of the target. In the
% frequency domain, the spectrum of the transmission signal is multiplied
% with the frequency response of the target.

%  filter your Tx signal with this frequency response for a given bubble ?--radius

%% Insert custom freq. response here!   
        end
    end
end

% Apply damping: 1/r (linear), 1/r^2 (cylindrical), 1/r^3 (spherical)
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
% Plot --------------------------------------------------------------------
figure;
subplot(211);
plot(tSim, rx(:, 1));
grid on;
title('Received signal');
subplot(212);
semilogx(tLagInMeters, corr);
hold on;
semilogx(tLagInMeters(loc), pk, 'rx');
grid on;
title(['Crosscorrelation: Transmit- & receive signal']);

%% Beamforming: Calculate array manifold vector (AMV)
NFFT = 2^nextpow2(2 * nRxSeqLength); 
%NFFT = nSig;
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
F=repmat(f,[NBeams, 1,NRx]);
% Reshape to [beams, NRx, freqs]
F = permute(F, [1 3 2]);
% Multiply tau with freq. supports
tauTimesf = bsxfun(@times,tTau', F);
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

%% Plot results as PPI (plan position indicator) plot
figure;
% Log of correlated output
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
