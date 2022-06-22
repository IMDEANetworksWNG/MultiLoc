clear 
clc
close all

att = 5e-1;

csi_data = +0.5i*ones(6,6,100) + (randn(6,6,100) + 1i*randn(6,6,100))*sqrt(2)*att;

% figure, plot(angle(squeeze(csi_data(1,1,:))))
% mean(angle(csi_data),3)
csi_data = sum(csi_data,3)/100;
% figure, plot(angle(squeeze(csi_data(1,1))))



% number of antennas
N = 6;

% frequency
freq = 60.48e9;

% speed of light
c = 3e8;

% the wavelength
lambda = c/freq;

% distance between antennas
d = lambda*0.58;

% step for the angle
step_angle = 1;

% take the codebook
[cb_az, theta_az] = Grid_AoA(step_angle, N,d,lambda);
[cb_el, theta_el] = Grid_AoA(step_angle, N,d,lambda);

csi_data(1,1) = 0;
csi_data(6,6) = 0;
csi_data(1,6) = 0;
csi_data(6,1) = 0;
csi_data(5,3) = 0;

% apply mD-track
[Az_estimated, El_estimated, att] = mD_track_2D(csi_data, cb_az, cb_el);

Az_estimated = rad2deg(theta_az(Az_estimated))
El_estimated = rad2deg(theta_el(El_estimated))
att