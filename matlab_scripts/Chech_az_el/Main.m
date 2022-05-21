clear 
close all
clc

% create the channel

% number of signals
n_s = 2;

% AoA
AoA = [-45,45];

% attenuation
alpha = [1 1];

% number of antennas
N = 6;

% frequency
freq = 60e9;

% speed of light
c = 3e8;

% the wavelength
lambda = c/freq;

% distance between antennas
d = lambda*0.5;

% number of samples;
n_samples = 1e3;
channel = zeros(N,n_samples);
for ii = 1:n_samples
    channel(:,ii) = channel_generator(alpha,AoA,N,d,lambda);
end

% step for the angle
step_angle = 1;

% create the codebo3ok
[cb_aoa, theta] = Grid_AoA(step_angle, N,d,lambda);

% create the correlation matrix
C = channel * channel';

% apply MUSIC
[ps_db, D] = MUSIC(C, cb_aoa, 1, 0);

figure, plot(rad2deg(theta),ps_db)