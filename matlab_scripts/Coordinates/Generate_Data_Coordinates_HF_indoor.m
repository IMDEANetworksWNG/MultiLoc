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
    angles_dp(ii) = possible_angles(index_min);
end
angles_dp = angles_dp*(-1);

angles_dp(:,1,:) = angles_dp(:,1,:) - aoa_offset(1);
angles_dp(:,2,:) = angles_dp(:,2,:) - aoa_offset(2);
angles_dp(:,3,:) = angles_dp(:,3,:) - aoa_offset(3);
angles_dp(:,4,:) = angles_dp(:,4,:) - aoa_offset(4);
angles_dp(:,5,:) = angles_dp(:,5,:) - aoa_offset(5);

angles_dp = angles_dp*(-1);
save("mat_files/indoor/HF/calculated_aoa_amb", "angles_dp")

%% Estimate the positions
estimated_point_x_hf = zeros(size(calculated_aoa));
estimated_point_y_hf = zeros(size(calculated_aoa));

for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points

            estimated_point_x_aux = sind(angles_dp(id_point,id_router, id_rotation))*estimated_distance(id_point,id_router, id_rotation);
            estimated_point_y_aux = cosd(angles_dp(id_point,id_router, id_rotation))*estimated_distance(id_point,id_router, id_rotation);


    %         figure, plot(estimated_point_x_p(id_point,id_router), estimated_point_y_p(id_point,id_router), "*")


            if (id_router == 1 || id_router == 3 || id_router == 4)
                estimated_point_x_hf(id_point,id_router, id_rotation) = RX_x(id_router) - (estimated_point_x_aux);
                estimated_point_y_hf(id_point,id_router, id_rotation) = RX_y(id_router) - (estimated_point_y_aux*(-1));
            else
                estimated_point_x_hf(id_point,id_router, id_rotation) = RX_x(id_router) - (estimated_point_x_aux*(-1));
                estimated_point_y_hf(id_point,id_router, id_rotation) = RX_y(id_router) - estimated_point_y_aux;
            end
        end
    end
end


save("mat_files/Coordinates/data_coordinates_hf_indoor","estimated_point_x_hf", "estimated_point_y_hf");

cd(pwd_str)