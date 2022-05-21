clear all
clc
close all


%% General stuff
pwd_str = pwd;

cd ../../../

addpath("functions")

load("mat_files/indoor_map/data_distance_angle_true.mat")
load("mat_files/indoor_map/points_coordinates.mat")
% load("mat_files/walls.mat")
load("mat_files/indoor_map/RX_coordinates.mat")
mkdir("mat_files/Coordinates/")


% RX_x = [RX_x_4, RX_x_5, RX_x_10];
% RX_y = [RX_y_4, RX_y_5, RX_y_10];



% angles_dp = angles_dp

% fig_map = Plot_Map_No_Points();


% load the ranging
load("mat_files/FTM/LF/estimated_distance_indoor.mat")
% load('mat_files/imdea_map/data_distance_angle_true.mat');

% estimated_distance = calculated_distance;

% load the angle
load("mat_files/mD_track/LF/Data_3D_Direct_Path_indoor.mat")

% angles_dp = angles_RX*(-1);
% estimated_distance = distances_RX;

% angles_dp(:,3) = angles_dp(:,3) - 6;
% angles_dp(:,2) = angles_dp(:,2) + 6;
% angles_dp([19,20,27:31],1) = angles_dp([19,20,27:31],1)*(-1);
[n_points, routers] = size(angles_dp);


% load the coordinates
load("mat_files/Coordinates/data_coordinates_lf_indoor");

error = sqrt((points_x - estimated_point_x_lf).^2 + (points_y - estimated_point_y_lf).^2);
error = error(:);

angles_dp = abs(angles_dp(:));

[angle_sort, index_sort] = sort(angles_dp);

error_sort = error(index_sort);

figure, plot(angle_sort, error_sort);

angle_th = 45;

error_below = error_sort(angle_sort < angle_th);
error_above = error_sort(angle_sort >= angle_th);

figure, cdfplot(error_below)
hold on
cdfplot(error_above)