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
time_end = 1;
bubbleTime = 3; % can be assumed as 0.1s
bubbleVelocity = 1.6/bubbleTime; % v = 1m / 0.1s;
for time_step = 1%:time_end 
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
tSig = 0.01; % in seconds
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

% fig=figure;
% subplot(211);
% grid on;
% plot(t, tx(:,1))
% title("Time domain signal") 
% best_plot_ever(fig)
% subplot(212);
% logTx = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
% grid on;
% plot(-f(1:NBins), logTx)
% title("Log. frequency spectrum ")
% best_plot_ever(fig)

%% Bubble environment settings
%  source location constrains a, b
x_lims=[-0.1 0.1];
y_lims=[20 20.2];
% y_lims=[1.2 1.4];
z_lims=[0 1];
Nbubbles=50;
bubbleOsc_lims = [-0.1,0.1];
% minRadius = 1000e-6;
% minAllowableDistance = max([585e-6, 2 * maxRadius]);
posTar = set_bubble_flare(x_lims, y_lims, z_lims, Nbubbles, bubbleVelocity, time_end, bubbleOsc_lims);
% posTar2 = set_bubble_flare([50 51], [30 32], [10 11], Nbubbles, bubbleVelocity, time_end, bubbleOsc_lims);
% posTar = [posTar1; posTar2];
% if time_step == 1
%     % Generate bubbles in some constrained space
%     posTarNew = set_bubble_source(x_lims, y_lims, z_lims, Nbubbles);
%     posTar = posTarNew;
% else
%     rng(time_step)
%     % posTarNew = set_bubble_source(x_lims, y_lims, z_lims, Nbubbles);
%     posTar(:,3) = posTar(:,3) + bubbleVelocity*ones(NTargets, 1);
%     bubbleOscillations = bubbleOsc_lims(1) + (bubbleOsc_lims(2) - bubbleOsc_lims(1))*rand(NTargets,2);
%     posTar(:,1:2) = posTar(:,1:2) + bubbleOscillations;
%     posTar = [posTar; posTarNew];
% end
NTargets = size(posTar, 1);
bDirectSound = 0;
%% Radius of bubbles
% Generate values from a normal distribution with mean and standard deviation
% a_range = linspace(8e-6,1000e-6,NTargets);
% a_range = linspace(585e-6,1000e-6,NTargets);
% radius_mean = 10e-3;
radius_std = 0.08e-3;%3e-3;
str_radius = num2str(radius_std);
str_radius = strrep(str_radius, '.', ',');
radius_mean = 8e-5;   % r= 8e-5 resonances in the 30-70 kHz
% radius_std = 0.8e-5;
a_range = radius_mean + radius_std.*randn(1,NTargets);                    %stochastic distribution
% a_range = radius_std*ones(1,NTargets);

%% Plot bubbles in 2D
x = posTar(:,1);
y = posTar(:,2);
z = posTar(:,3);
bubbles_mov = figure;
% hold on
plot3(x,y,z, '-ok',"MarkerEdgeColor" ,	"#4DBEEE")
axis_x = x_lims+bubbleOsc_lims*time_end*bubbleVelocity;
axis_y = y_lims+bubbleOsc_lims*time_end*bubbleVelocity;
axis([axis_x   axis_y     0 time_end*bubbleVelocity])
axis equal
grid on
view([0  90])
% refreshdata
% drawnow
% Frame = getframe(bubbles_mov);
% make_gif(Frame, time_step, "Bubble_mov2.gif");
%% Draw circles of different radii
% n=NTargets;
% R=radius_std*ones(n,1)*40;
% P=[0:0.1:2*pi 0];
% xr=R.*cos(P);
% yr=R.*sin(P);
% X_C=bsxfun(@plus,x,xr);
% Y_C=bsxfun(@plus,y,yr);
% fig = figure;
% plot(X_C',Y_C','b')
% xlabel('X');ylabel('Y')
% axis equal
% saveas(gca, "thesis_pics/bubbles_plot_2D","png")
%% SONAR-sytem element positions (line arrays)
% Uniform line array element equidistant element spacing 
% around array center (spaced on x-axis)
rxCenterElement = (NRx + 1) / 2;
txCenterElement = (NTx + 1) / 2;
turnRx = [cosd(0) sind(0) sind(0)];
turnTx = [cosd(0) sind(0) sind(0)];
posRx = (((1:NRx) - rxCenterElement) .* dRx)' .*turnRx+ centerRx;
posTx = (((1:NTx) - txCenterElement) .* dTx)' .*turnTx+ centerTx;
% figure;plot(posRx(:,1:2))
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

%% Calculate propagation delays between bubbles
tBubblePropagationDist = zeros(NTargets, NTargets);
% Propagation distance: Tx -> target -> Rx
for iTx = 1:NTx
for iTar1 = 1:NTargets
    for iTar2 = iTar1+1:NTargets
        dTar1ToTar2 = norm(posTar(iTar1, :) - posTar(iTar2, :));
        tBubblePropagationDist(iTar1, iTar2) = dTar1ToTar2;
    end
end
end
tBubblePropagationDist=tBubblePropagationDist'+triu(tBubblePropagationDist',1)';
% Propagation distance -> Propagation time
tBubblePropagationTime = tBubblePropagationDist ./ cWater; % in seconds
tBubblePropagationTime = tBubblePropagationTime .* fs; % in samples

%% Calculate received signals
% Geometric spreading loss damping
geoSpreadLoss = 2;
% Max rx. sequence length (signal duration + max propagation time)
nRxSeqLength = nSig + ceil(max(tPropagationTime(:))) + ceil(max(tBubblePropagationTime(:)));
% nRxSeqLength_bub = nRxSeqLength + ceil(max(tBubblePropagationTime(:)));
rx = zeros(nRxSeqLength, NRx);
rx_multi = zeros(nRxSeqLength, NRx);

%% Background noise
noise_level_dB = -65;
noise_level_linear = 10^(noise_level_dB/10);
noise_add = randn(nRxSeqLength, NRx) * noise_level_linear; 
% rx = rx + noise_add; rx_multi = rx_multi + noise_add;

%% Bubbles freq. response 
sigma_bs = bubble_response_model(f,a_range,1);
f_sigma_bs = abs(Tx(:, iTx)).*sigma_bs(1:NBins,:);
theta = angle(Tx(:, iTx));
mixedResponsesFreq_init = f_sigma_bs.*exp(i*theta);
mixedResponsesTime_init = ifft(mixedResponsesFreq_init, NFFT,'symmetric');
mixedResponsesTime_init = mixedResponsesTime_init(1:nSig,:);
%% 
% Plot --------------------------------------------------------------------
% fig = figure; 
% subplot(211);
% for sb_i = 1: NTargets
%     hold on 
%     logBs = 10*log10(sigma_bs(1:NBins,sb_i));
%     plot(-f(1:NBins), logBs)
%     title("Log. frequency spectrum, only bubbles")
%     best_plot_ever(fig)
%     ylim([-200 0])
%     hold off
% end
% subplot(212);
% % fig = figure; 
% logFs = 20*log10(abs(f_sigma_bs(:,:))./max(abs(f_sigma_bs(:,:))));
% plot(-f(1:NBins), logFs(:,10))
% title("Log. frequency spectrum, with a bubble")
% ylim([-200 0])
% best_plot_ever(fig)
%% Add bubble response to the received signal
for iTx = 1:NTx
    for iRx = 1:NRx
        for iTar = 1:NTargets + bDirectSound
            iStart = floor(tPropagationTime(iTx, iTar, iRx))+1;
            iEnd = floor(iStart + length(tx(:, iTx))) - 1;
%% Custom freq. response here! 
            % Add reflected tx-signal at delay time to rx_-
            % rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) +  tx(:, iTx);
            tarList = 1:NTargets;
            tarList(iTar) = [];
            for iTar2 = iTar+1:NTargets
                iStartBubble = floor(tPropagationTime(iTx, iTar, iRx)+tBubblePropagationTime(iTar, iTar2))+1;
                iEndBubble = floor(iStart + length(tx(:, iTx))) - 1;
                mixed_resp = f_sigma_bs(:,iTar).*mixedResponsesFreq_init(:,iTar2);
                mixed_resp = ifft(mixed_resp, NFFT,'symmetric');
                mixed_resp=mixed_resp(1:nSig,:);
                rx_multi(iStartBubble:iEndBubble, iRx) = rx_multi(iStartBubble:iEndBubble, iRx) + mixed_resp;
            end
%             idea: to add phase shift (imag part of the Tx to the mixed
%             freq resp)
            rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) + mixedResponsesTime_init(:,iTar);
            rx_multi(iStart:iEnd, iRx) = rx_multi(iStart:iEnd, iRx) + mixedResponsesTime_init(:,iTar);
            % rx(iStart:iEnd, iRx) = rx(iStart:iEnd, iRx) + mixed_resp(1:nSig, iTx);

  
        end
    end
end
%% Reconstruction of the signal
Y = fft(rx, NFFT);
Y = Y(1:NBins, :);
Y_multi = fft(rx_multi, NFFT);
Y_multi = Y_multi(1:NBins, :);

f_half =  -f(1:NBins);
[~, ffmin] = min(abs(f_half-fMin));
[~, ffmax] = min(abs(f_half-fMax));

H_hat = abs(Y(:,10)) ./ abs(Tx(:,1));
H_hat = H_hat.*hamming(length(H_hat));
H_multi_hat = abs(Y_multi(:,10)) ./ abs(Tx(:,1));
H_multi_hat = H_multi_hat.*hamming(length(H_multi_hat));

% H_hat(ffmin:ffmax) = abs(Y(ffmin:ffmax,10)) ./ abs(Tx(ffmin:ffmax,1));
% H_multi_hat(ffmin:ffmax) = abs(Y_multi(ffmin:ffmax,10)) ./ abs(Tx(ffmin:ffmax,1));

% 
% fig=figure;
% % subplot(211)
% plot(-f(1:NBins), 10*log10(abs(H_hat(:,1))./max(abs(H_hat(ffmax:ffmin,1)))),'LineWidth', 1.5);
% hold on
% plot(-f(1:NBins), 10*log10(abs(H_multi_hat(:,1))./max(abs(H_multi_hat(ffmax:ffmin,1)))));
% legend("$\hat{H}$", "$\hat{H}_{\mathrm{multi}}$")
% xlabel('Frequency (Hz)');ylabel('Magnitude, dB'); 
% % xlim([f_half(ffmin) f_half(ffmax)])
% best_plot_ever(fig)
% saveas(gca, "thesis_pics/multi_scat_50kHz_r"+str_radius+"_Hhat_norm_comparison","png");

fig=figure;
% subplot(211)
plot(-f(1:NBins),abs(H_hat),'LineWidth', 1.5);
hold on
plot(-f(1:NBins), abs(H_multi_hat),'-');
legend("$\hat{H}$", "$\hat{H}_{\mathrm{multi}}$")
xlabel('Frequency (Hz)'); ylabel('Amplitude'); 
% xlim([f_half(ffmin) f_half(ffmax)])
best_plot_ever(fig)
% saveas(gca, "thesis_pics/multi_scat_50kHz_r"+str_radius+"_Hhat_amplitude_comparison","png");
%% Interval amplitude plots H_hat

% (ffmin:ffmax) 
% fig=figure;
% % subplot(211)
% % plot(-f(ffmin:ffmax) , 20*log10(abs(H_hat(:,1))./max(abs(H_hat(ffmin:ffmax) ))),'LineWidth', 1.5);
% hold on
% % ranged_H_hat = 20*log10(abs(H_multi_hat(:,1))./max(abs(H_multi_hat(ffmin:ffmax) )));
% % plot(-f(1:NBins), ranged_H_hat);
% legend("$H_{hat}$", "$H_{multi hat}$")
% xlabel('Frequency (Hz)');
% best_plot_ever(fig)
% saveas(gca, "thesis_pics/multi_scat_50kHz_r"+str_radius+"_Hhat_norm_comparison","png");


%% Apply damping: 1/r (linear), 1/r^2 (cylindrical), 1/r^3 (spherical)
propagationTimesPerSample = (0:nRxSeqLength)./ fs;
propagationDistancePerSample = propagationTimesPerSample * cWater;
damping = 1./(propagationDistancePerSample.^geoSpreadLoss);
% damping1 = repmat(damping, NRx,1);
% rx = rx - damping1';
% rx_multi = rx_multi - damping1';
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
%% Plot --------------------------------------------------------------------
fig=figure;
subplot(211);
plot(tSim, rx(:, 1));
grid on;
title('Received signal');
best_plot_ever(fig)
subplot(212);
semilogx(tLagInMeters, corr);
hold on;
semilogx(tLagInMeters(loc), pk, 'rx');
grid on;
title('Crosscorrelation: Transmit- and receive signal');
best_plot_ever(fig)
% saveas(gca, "thesis_pics/multi_scat_50kHz_r"+str_radius+"_corr-distance","png");


%% Beamforming: Calculate array manifold vector (AMV)
NFFT = 2^nextpow2(2 * nRxSeqLength); 
% NFFT = nSig;
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

% rx_multi signal to freq. domain
Rx_multi = fft(rx_multi, NFFT);
Rx_multi = Rx_multi(1:NBins, :);

% tx signal to freq. domain
Tx = fft(tx, NFFT);
Tx = Tx(1:NBins, :);
% Matched filtering: Find tx-sequence in each rx-sequence
Mf_multi = conj(Tx) .* Rx_multi;

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
H_hat = abs(Y(:,round(NRx/2))) ./ abs(Tx(:,1));
% H_hat = H_hat.*hamming(length(H_hat));
% Normalize to max
logH_hat = 20*log10(abs(H_hat(:,1))./max(abs(H_hat(:,1))));

Y_multi = Mf_multi;
H_hat_multi = abs(Y_multi(:,round(NRx/2))) ./ abs(Tx(:,1));
% H_hat_multi = H_hat_multi.*hamming(length(H_hat_multi));
% Normalize to max
logH_hat_multi = 20*log10(abs(H_hat_multi(:,1))./max(abs(H_hat_multi(:,1))));

f = linspace(-fs/2, fs/2, NFFT);
f_half =  -f(1:NBins);

logBs = 10*log10(sigma_bs(1:NBins,1));

f = linspace(-fs/2, fs/2, NFFT);
fig=figure;
subplot(311)
% logRx_noB = X(:,1);
logTx_noB = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
plot(-f(1:NBins), logTx_noB(:, 1));
ylim([-100 0])
grid on;
title('Transmitted, spectrum');
best_plot_ever(fig)

subplot(312)
% logRx_withB = Y(:,1);
logRx_withB = 20*log10(abs(Y(:,1))./max(abs(Y(:,1))));
logRx_multi_withB = 20*log10(abs(Y_multi(:,1))./max(abs(Y_multi(:,1))));
plot(-f(1:NBins), logRx_withB(:, 1));
ylim([-100 0])
grid on;
hold on
% plot(-f(1:NBins), logRx_multi_withB(:, 1));
title('Received correlated, spectrum');
ylabel('Magnitude, dB')
best_plot_ever(fig)

subplot(313)

hold on
% plot(2*pi/c*-f(1:NBins)*radius_std, logH_hat(:, 1));
plot(-f(1:NBins), logH_hat);
% plot(-f(1:NBins), logH_hat_multi);
plot(-f(1:NBins), logBs, '-.','LineWidth',1.5);
ylim([-100 0])
legend("reconstruction", "bubble", 'Location','best')
title('Reconstructed correlated, spectrum');
xlabel('Frequency (Hz)');
best_plot_ever(fig)
% saveas(gca, "thesis_pics/multi_scat_50kHz_r"+str_radius+"_reconstruct","png");
% saveas(gca, "thesis_pics/multi_scat_reconstr_corr-set_50kHz_r"+num2str(radius_b),"png");


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
% axis_limit = 110;
% clim([-40 0])
% xlim([-axis_limit axis_limit]); 
% ylim([0 axis_limit])
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
% saveas(gca, "thesis_pics/multi_scat_50kHz_r"+str_radius+"_sonar-plot","png");
end

%% Functions
function posTar = set_bubble_source(x_lims, y_lims, z_lims, Nbubbles)
    posTarX = x_lims(1) + (x_lims(2)-x_lims(1))*rand(Nbubbles,1);
    posTarY = y_lims(1) + (y_lims(2)-y_lims(1))*rand(Nbubbles,1);
    posTarZ = z_lims(1) + (z_lims(2)-z_lims(1))*rand(Nbubbles,1);
    posTar = [posTarX posTarY posTarZ];
end
function posTar = set_bubble_flare(x_lims, y_lims, z_lims, Nbubbles, bubbleVelocity, time_end, bubbleOsc_lims)
    for tt =1:time_end
        rng(tt)
        if tt == 1
            posTarNew = set_bubble_source(x_lims, y_lims, z_lims, Nbubbles);
            posTar = posTarNew;
        else
            NTargets = size(posTar, 1);
            posTar(:,3) = posTar(:,3) + bubbleVelocity*ones(NTargets, 1);
            bubbleOscillations = bubbleOsc_lims(1) + (bubbleOsc_lims(2) - bubbleOsc_lims(1))*rand(NTargets,2);
            posTar(:,1:2) = posTar(:,1:2) + bubbleOscillations;
            posTar = [posTar; posTarNew];
        end
    end
end
function make_gif(Frame, ii, filename)
    im = frame2im(Frame); 
    [imind, CM] = rgb2ind(im,256); 
    % Write the animation to the gif File: MYGIF 
    if ii == 1 
      imwrite(imind, CM,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind, CM,filename,'gif','WriteMode','append'); 
    end 
end

