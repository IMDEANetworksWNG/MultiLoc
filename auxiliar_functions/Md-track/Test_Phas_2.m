clear
clc
close all

load("pan_0_tilt_0.mat")

% move the phase to radians
phases = (phases/1024)*2*pi;

% move to complex
csi_data = exp(1i*phases);


% take the data for antenna 2
csi_data_test = csi_data(:,17);
figure, plot(angle(csi_data_test))

[csi_data_clean,offset] = Sanitize(csi_data_test);
figure, plot(angle(csi_data_clean))
% figure, plot(angle(csi_data_clean*exp(1i*offset)))