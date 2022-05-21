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
load("mat_files/FTM/LF/estimated_distance.mat")

% load the angle
load("mat_files/mD_track/LF/Data_3D_Direct_Path.mat")

angles_dp(1:74,2) = angles_dp(1:74,2) - 3.5;
angles_dp(26:38,5) = angles_dp(26:38,5) - 5;
angles_dp(1:74,3) = angles_dp(1:74,3) - 7;
% angles_dp(92,2) = 0;
% angles_dp(93,2) = 8;
% angles_dp(77,5) = 0;
% 
% estimated_distance(63,1) = 6.5;
% estimated_distance(106,1) = 16.5;


%% Estimate the positions
for id_router = 1:routers
    for id_point = 1:n_points

        estimated_point_x_p(id_point,id_router) = sind(angles_dp(id_point,id_router))*estimated_distance(id_point,id_router);
        estimated_point_y_p(id_point,id_router) = cosd(angles_dp(id_point,id_router))*estimated_distance(id_point,id_router);


%         figure, plot(estimated_point_x_p(id_point,id_router), estimated_point_y_p(id_point,id_router), "*")


        if (id_router == 1 || id_router == 3 || id_router == 4)
            estimated_point_x_lf(id_point,id_router) = RX_x(id_router) - (estimated_point_x_p(id_point,id_router));
            estimated_point_y_lf(id_point,id_router) = RX_y(id_router) - (estimated_point_y_p(id_point,id_router)*(-1));
        else
            estimated_point_x_lf(id_point,id_router) = RX_x(id_router) - (estimated_point_x_p(id_point,id_router)*(-1));
            estimated_point_y_lf(id_point,id_router) = RX_y(id_router) - estimated_point_y_p(id_point,id_router);
        end
    end
end



save("mat_files/Coordinates/data_coordinates_lf_old","estimated_point_x_lf", "estimated_point_y_lf");

cd(pwd_str)