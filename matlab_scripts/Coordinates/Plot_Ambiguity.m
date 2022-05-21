clear 
close all
clc

% create the grid for the angles
angles = -90:0.1:90;
% calculate the phase
phases = 2*0.58*pi*sin(deg2rad(angles));


% plot the phases for the angles
figure, plot(angles,phases)
hold on
plot(angles,ones(length(angles),1)*pi)
plot(angles,ones(length(angles),1)*-pi)

% center everything between +/- pi
phases(phases > pi) = phases(phases > pi) - 2*pi;
phases(phases < -pi) = phases(phases < -pi) + 2*pi;

% plot it again
figure, plot(angles,phases)
hold on
plot(angles,ones(length(angles),1)*pi)
plot(angles,ones(length(angles),1)*-pi)