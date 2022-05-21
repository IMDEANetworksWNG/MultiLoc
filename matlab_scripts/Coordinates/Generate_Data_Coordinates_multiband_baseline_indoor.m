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


[n_points, routers, rotations] = size(estimated_distance);

% load the angle
load("mat_files/Jade/Data_MUSIC_AoA_Direct_Path_3.mat")
angles_dp_lf = angles_dp;
angles_dp_lf(1:74,2) = angles_dp_lf(1:74,2) - 3.5;
angles_dp_lf(26:38,5) = angles_dp_lf(26:38,5) - 5;
angles_dp_lf(1:74,3) = angles_dp_lf(1:74,3) - 7;


%% Estimate the positions
estimated_point_x_mb_bs = zeros(size(estimated_distance));
estimated_point_y_mb_bs = zeros(size(estimated_distance));

for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points

            estimated_point_x_aux = sind(angles_dp(id_point,id_router))*estimated_distance(id_point,id_router, id_rotation);
            estimated_point_y_aux = cosd(angles_dp(id_point,id_router))*estimated_distance(id_point,id_router, id_rotation);


    %         figure, plot(estimated_point_x_p(id_point,id_router), estimated_point_y_p(id_point,id_router), "*")


            if (id_router == 1 || id_router == 3 || id_router == 4)
                estimated_point_x_mb_bs(id_point,id_router, id_rotation) = RX_x(id_router) - (estimated_point_x_aux);
                estimated_point_y_mb_bs(id_point,id_router, id_rotation) = RX_y(id_router) - (estimated_point_y_aux*(-1));
            else
                estimated_point_x_mb_bs(id_point,id_router, id_rotation) = RX_x(id_router) - (estimated_point_x_aux*(-1));
                estimated_point_y_mb_bs(id_point,id_router, id_rotation) = RX_y(id_router) - estimated_point_y_aux;
            end
        end
    end
end


save("mat_files/Coordinates/data_coordinates_multiband_baseline_indoor","estimated_point_x_mb_bs", "estimated_point_y_mb_bs");

cd(pwd_str)