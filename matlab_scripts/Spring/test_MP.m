clear 
% close all
% clc

% create the channel in time
K = 256;
channel = zeros(K,1);
channel(3) = (1 + 1i) /sqrt(2);

% channel in freq
channel_freq = fft(channel);
channel_freq = channel_freq(6:(256-6));

load("csi_toa_interp.mat")

channel_freq = csi_toa_interp;
% the pencil
P = 100;
[z_l] = matrix_pencil(channel_freq,P);
estimated_delay = (K/(2*pi)).*(angle((z_l)))*12.5

[S_toa, times] = ToA_Phases(0.1, K, 80, 1:245);
spectrum_u = S_toa' * channel_freq;
figure, plot(times,abs(spectrum_u))
