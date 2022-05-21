clear
close all
clc

pwd_str = pwd;

cd ../../

% load the true data
load("mat_files/indoor_map/data_distance_angle_true.mat")
load("mat_files/indoor_map/points_coordinates.mat")
load("mat_files/indoor_map/RX_coordinates.mat")

% load the map
openfig("mat_files/indoor_map/map_no_points.fig")
hold on
plot(RX_x(5), RX_y(5), "r*", "Linewidth", 4)
plot(points_x(1), points_y(1), "b*", "Linewidth", 4)

% load HF data
load("mat_files/indoor/HF/CSI/csi_indoor.mat")
load("mat_files/indoor/HF/CSI/aoa_optimal_rotation.mat")
load("mat_files/indoor/HF/FTM/ftm_indoor.mat")

rotation = 5;

angle_hf = squeeze(calculated_aoa(1,5,rotation)) + 4
distance_hf = squeeze(calculated_distance(1,5,rotation))
angle_hf = angle_hf*(-1);
% angle_hf = 0;

plot(RX_x(5)+[(0:0.1:15)*sind(angle_hf)], RX_y(5)+[(0:0.1:15)*cosd(angle_hf)*(-1)], '-', 'LineWidth', 2);
