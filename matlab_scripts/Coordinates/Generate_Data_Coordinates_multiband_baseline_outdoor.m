clear all
clc
%close all


%% General stuff
pwd_str = pwd;

cd ../../

addpath("functions")

load("mat_files/outdoor_map/data_distance_angle_true.mat")
load("mat_files/outdoor_map/points_coordinates.mat")
% load("mat_files/walls.mat")
load("mat_files/outdoor_map/RX_coordinates.mat")
mkdir("mat_files/Coordinates/")


% RX_x = [RX_x_4, RX_x_5, RX_x_10];
% RX_y = [RX_y_4, RX_y_5, RX_y_10];



% angles_dp = angles_dp

% fig_map = Plot_Map_No_Points();


% load the ranging
load("mat_files/outdoor/HF/FTM/ftm_outdoor.mat")
load('mat_files/outdoor_map/data_distance_angle_true.mat');

estimated_distance = calculated_distance;

% load the angle
load("mat_files/Jade/Data_MUSIC_AoA_Direct_Path_Outdoor_3.mat")

% Calibrate
angles_dp(:,2) = angles_dp(:,2) + 3.5;
angles_dp(:,3) = angles_dp(:,3) - 2;

load('mat_files/outdoor/HF/optimal_rotation_index.mat')

[n_points, routers, rotations] = size(calculated_distance);


%% Localize: Estimate the positions

estimated_direct_path_angle    = angles_dp;% - median(abs(calculated_aoa - repmat(angles_RX, [1, 1, 8])), 'omitnan');
estimated_direct_path_distance = calculated_distance;

estimated_direct_path_angle(:, [1 4])    = estimated_direct_path_angle(:, [4 1]);


% https://en.wikipedia.org/wiki/Rotation_of_axes
for point_id=1:16
    for ap_id=1:4
        for rotation_id=1:8
            
            estimated_point_x_p(point_id, ap_id, rotation_id) = sind(estimated_direct_path_angle(point_id,ap_id))*estimated_direct_path_distance(point_id,ap_id, rotation_id);
            estimated_point_y_p(point_id, ap_id, rotation_id) = cosd(estimated_direct_path_angle(point_id,ap_id))*estimated_direct_path_distance(point_id,ap_id, rotation_id);

            % Rotate respect the AP
            x_center = RX_x(ap_id);
            y_center = RX_y(ap_id);

            if (ap_id == 1) % 43

                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x_mb_bs(point_id,ap_id, rotation_id)));
                theta = deg2rad(225);
                v = [estimated_point_x_mb_bs(point_id,ap_id, rotation_id);estimated_point_y_mb_bs(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 2) % 44

                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x_mb_bs(point_id,ap_id, rotation_id)));
                theta = deg2rad(-45);
                v = [estimated_point_x_mb_bs(point_id,ap_id, rotation_id);estimated_point_y_mb_bs(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;   
                vo = so + center; 
                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 3) %45

                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x_mb_bs(point_id,ap_id, rotation_id)));
                theta = deg2rad(45);
                v = [estimated_point_x_mb_bs(point_id,ap_id, rotation_id);estimated_point_y_mb_bs(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 4) % 46

                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x_mb_bs(point_id,ap_id, rotation_id)));
                theta = deg2rad(-90);
                v = [estimated_point_x_mb_bs(point_id,ap_id, rotation_id);estimated_point_y_mb_bs(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x_mb_bs(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y_mb_bs(point_id,ap_id, rotation_id) = vo(2,:);
            end
            
%             estimated_point_x_mb_bs(point_id, ap_id, rotation_id);
%             estimated_point_y_mb_bs(point_id, ap_id, rotation_id);
              plot(estimated_point_x_mb_bs(point_id, ap_id, rotation_id), estimated_point_y_mb_bs(point_id, ap_id, rotation_id), "*")

        end
    end
end


save("mat_files/Coordinates/data_coordinates_multiband_baseline_outdoor","estimated_point_x_mb_bs", "estimated_point_y_mb_bs");

cd(pwd_str)