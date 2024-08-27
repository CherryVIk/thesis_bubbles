% models_trial
% Thuraisingham 1
% Anderson 2
% Church 3

clc
clear
close all

f_range = linspace(0.1,300,30000)*1000;
% fs = 192000;
% NFFT = 32768;
% NBins = NFFT / 2 + 1;
% f_range = linspace(-fs/2, fs/2, NFFT);%NFFT
% a_range = linspace(8e-4,4e-3,10); % very slow for high number of radiuses
% at Anderson model
a_range = 3e-3;% 0.3e-3; %a_range(aa); % bubble radius (m)
TS_thur = 10*log10(bubble_response_model(f_range,a_range, 1));
TS_and = 10*log10(bubble_response_model(f_range,a_range, 2));
% TS_church = 10*log10(bubble_response_model(f_range,a_range, 3));
% Plot log(ka) x TS
a = a_range(1);
kk = 1; % at specific radius
c=1500; % the speed of sound in water
ka = 2*pi/c*f_range'*a_range;
fig=figure;

semilogx(ka(:,kk), TS_thur(:,kk),'LineWidth', 1.5, 'DisplayName','Thuraisingham')
hold on
semilogx(ka(:,kk), TS_and(:,kk),'LineWidth', 1.5,'LineStyle', '-.','DisplayName','Anderson')
ylim([-100 0])
legend('Thuraisingham','Anderson')
xlabel('log(ka)');ylabel('TS (dB re 1 $m^2$)')
titlename = "Plot ka x TS. TS for a sphere with a=" + (a*100) + " cm";
title(titlename)
best_plot_ever(fig)
% saveas(gca, "thesis_pics/plot_thur_ander_50kHz_a=" + num2str(a*100) + "cm","png");
%% Plot freq x TS
a = a_range(1);
fig=figure;
% subplot(311)
hold on
plot(f_range/1000, TS_thur);
plot(f_range/1000, TS_and);
xlabel('Freq (kHz)');ylabel('TS (dB re 1 m^2)')
titlename = "Plot freq x TS. TS for a sphere with a=" + num2str(a*100) + " cm";
title(titlename)
subtitlename = "Thuraisingham";
subtitle(subtitlename)
best_plot_ever(fig)
%% 

% Plot freq x TS
a = a_range(1);
fig=figure;
% subplot(311)
hold on
plot(-f_range(1:end/2+1)/1000, TS_and(1:NBins),'DisplayName','Anderson');
plot(-f_range(1:end/2+1)/1000, TS_thur(1:NBins),'DisplayName','Thuraisingham');
xlabel('Freq (kHz)');ylabel('TS (dB re 1 m^2)')
titlename = "Plot freq(end/2:end) x TS. TS(1:NBins) for a sphere with a=" + num2str(a*100) + " cm";
title(titlename)

% subtitlename = "Thuraisingham";
% subtitle(subtitlename)
best_plot_ever(fig)
% best_plot_ever(fig)
% saveas(gca, "thesis_pics/_50kHz_a=" + num2str(a*100) + "cm","png");
%% 
subplot(312)
hold on
plot(f_range/1000, TS_and);
xlabel('Freq (kHz)');ylabel('TS (dB re 1 m^2)')
subtitlename = "Anderson";
subtitle(subtitlename)
best_plot_ever(fig)

subplot(313)
hold on
plot(f_range/1000, TS_church);
xlabel('Freq (kHz)');ylabel('TS (dB re 1 m^2)')
subtitlename = "Church";
subtitle(subtitlename)
%% Plot ka x TS (three plots)
c=1500;
ka = 2*pi/c*f_range'*a_range;
kk = 1; % at specific radius
figure;
% subplot(211)

semilogx(ka(:,kk), TS_thur(:,kk))
hold on
semilogx(ka(:,kk), TS_and(:,kk))
semilogx(ka(:,kk), TS_church(:,kk))
ylim([-100 0])
legend('Thuraisingham','Anderson','Church')
xlabel('log(ka)');ylabel('TS (dB re 1 m^2)')
titlename = "Plot ka x TS. TS for a sphere with a=" + (a*100) + " cm";
title(titlename)

%% Functions
function best_plot_ever(fig)
    font = 20;
    a = 36; % set this parameter and keep it forever
    b = 0.55; % feel free to play with this ratio
    set(gca,'FontSize',font)
    set(findall(fig,'-property','Box'),'Box','off') % optional
    set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
    set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(fig,'Units','centimeters','Position',[3 3 a b*a])
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    set(gca, 'XDir', 'normal', 'YDir', 'normal');
    grid on
end