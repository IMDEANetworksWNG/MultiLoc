clear all
clc
% close all


%% General stuff
pwd_str = pwd;

cd ../../

addpath("functions")

load("mat_files/indoor_map/data_distance_angle_true.mat")
load("mat_files/indoor_map/points_coordinates.mat")
% load("mat_files/walls.mat")
load("mat_files/indoor_map/RX_coordinates.mat")
mkdir("mat_files/Coordinates/")


% RX_x = [RX_x_4, RX_x_5, RX_x_10];
% RX_y = [RX_y_4, RX_y_5, RX_y_10];

routers_csi = [4, 5,6,7, 10];



routers = 5;
n_points = 110;

% angles_dp = angles_dp

% fig_map = Plot_Map_No_Points();


% load the ranging
load("mat_files/indoor/HF/FTM/ftm_indoor.mat")
estimated_distance = calculated_distance;
load("mat_files/indoor/HF/FTM/ftm_optimal_rotation")
load("mat_files/indoor/HF/CSI/aoa_optimal_rotation.mat")

% load the angle
load("mat_files/indoor/HF/CSI/csi_indoor.mat")
% load the offset
load("mat_files/indoor/HF/CSI/aoa_offset.mat")
% angles_dp = repmat(angles_RX,1,1,8)*(-1);





[n_points, routers, rotations] = size(calculated_aoa);

% load the angle
load("mat_files/mD_track/LF/Data_3D_Direct_Path.mat")
angles_dp_lf = angles_dp;
angles_dp_lf(1:74,2) = angles_dp_lf(1:74,2) - 3.5;
angles_dp_lf(26:38,5) = angles_dp_lf(26:38,5) - 5;
angles_dp_lf(1:74,3) = angles_dp_lf(1:74,3) - 7;


%% remove ambiguity
angles_dp = calculated_aoa;
angles_dp = angles_dp*(-1);

% load the index to the optimal 
load("mat_files/indoor/HF/optimal_rotation_index.mat")

angles_optimal = zeros(110,5);
for id_point = 1:n_points
    for id_router = 1:routers
        index_optimal = indoor_optimal_rotation_index(id_point,id_router);
        angles_optimal(id_point, id_router) = angles_dp(id_point, id_router,index_optimal);
        
    end
end

angles_dp_2 = angles_dp;

max_angle = abs(rad2deg(asin((0.58*2*pi - 2*pi)/(0.58*2*pi))));

index_amb = find((abs(angles_dp) - max_angle) > 0);

for ii = index_amb.'
    angles_dp_aux = angles_dp(ii);
    
    phase_value = sind(angles_dp_aux)*0.58*2*pi;
    
    if(phase_value >= 0)
        phase_value = phase_value - 2*pi;

    else
        phase_value = phase_value + 2*pi;
    end
    
    angles_dp_aux_amb = rad2deg(asin(phase_value/(0.58*2*pi)));
    
    possible_angles = [angles_dp_aux angles_dp_aux_amb];
    
    [row, column, ~] = ind2sub(size(angles_dp), ii);
    
    angle_dp_lf_aux = angles_dp_lf(row, column);
    
    [~, index_min] = min(abs(angle_dp_lf_aux - possible_angles));
    angles_dp_2(ii) = possible_angles(index_min);
end

% angles_dp(:,1,:) = angles_dp(:,1,:) - aoa_offset(1);
% angles_dp(:,2,:) = angles_dp(:,2,:) - aoa_offset(2);
% angles_dp(:,3,:) = angles_dp(:,3,:) - aoa_offset(3);
% angles_dp(:,4,:) = angles_dp(:,4,:) - aoa_offset(4);
% angles_dp(:,5,:) = angles_dp(:,5,:) - aoa_offset(5);




angles_optimal_amb = zeros(110,5);
for id_point = 1:n_points
    for id_router = 1:routers
        index_optimal = indoor_optimal_rotation_index(id_point,id_router);
        angles_optimal_amb(id_point, id_router) = angles_dp(id_point, id_router,index_optimal);
        
    end
end

angles_optimal - angles_optimal_amb

cd(pwd_str)