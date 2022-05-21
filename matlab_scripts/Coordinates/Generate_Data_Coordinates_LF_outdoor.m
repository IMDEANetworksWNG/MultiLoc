clear all
clc
% close all


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

%routers_csi = [4, 5,6,7, 10];



% angles_dp = angles_dp

% fig_map = Plot_Map_No_Points();


% load the ranging
load('mat_files/FTM/LF/estimated_distance_outdoor.mat')


% load the angle
load("mat_files/mD_track/LF/Data_3D_Direct_Path_outdoor.mat")

% Calibrate
angles_dp(:,2) = angles_dp(:,2) + 3.5;
angles_dp(:,3) = angles_dp(:,3) - 2;
%% Localize: Estimate the positions

uiopen('mat_files/outdoor_map/outdoor_map.fig',1);
hold on
pause(1)

estimated_direct_path_angle    = angles_dp;
estimated_direct_path_distance = estimated_distance;

estimated_direct_path_angle(:, [1 4])    = estimated_direct_path_angle(:, [4 1]);
estimated_direct_path_distance(:, [1 4]) = estimated_direct_path_distance(:, [4 1]);
%RX_x(:, [1 4]) = RX_x(:, [4 1]);
%RX_y(:, [1 4]) = RX_y(:, [4 1]);

% https://en.wikipedia.org/wiki/Rotation_of_axes
for point_id=1:16
    for ap_id=1:4
            
        estimated_point_x_p(point_id, ap_id) = sind(estimated_direct_path_angle(point_id,ap_id))*estimated_direct_path_distance(point_id,ap_id);
        estimated_point_y_p(point_id, ap_id) = cosd(estimated_direct_path_angle(point_id,ap_id))*estimated_direct_path_distance(point_id,ap_id);

        % Rotate respect the AP
        x_center = RX_x(ap_id);
        y_center = RX_y(ap_id);

        if (ap_id == 1) % 43

            estimated_point_x_lf(point_id,ap_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id));
            estimated_point_y_lf(point_id,ap_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id)*(-1));

            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x_lf(point_id,ap_id)));
            theta = deg2rad(45);
            v = [estimated_point_x_lf(point_id,ap_id);estimated_point_y_lf(point_id,ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;
            vo = so + center; 
            estimated_point_x_lf(point_id,ap_id) = vo(1,:);
            estimated_point_y_lf(point_id,ap_id) = vo(2,:);

        elseif (ap_id == 2) % 44

            estimated_point_x_lf(point_id,ap_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id)*(-1));
            estimated_point_y_lf(point_id,ap_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id));

            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x_lf(point_id,ap_id)));
            theta = deg2rad(-45);
            v = [estimated_point_x_lf(point_id,ap_id);estimated_point_y_lf(point_id,ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;   
            vo = so + center; 
            estimated_point_x_lf(point_id,ap_id) = vo(1,:);
            estimated_point_y_lf(point_id,ap_id) = vo(2,:);

        elseif (ap_id == 3) %45

            estimated_point_x_lf(point_id,ap_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id)*(-1));
            estimated_point_y_lf(point_id,ap_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id));

            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x_lf(point_id,ap_id)));
            theta = deg2rad(35);
            v = [estimated_point_x_lf(point_id,ap_id);estimated_point_y_lf(point_id,ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;
            vo = so + center; 
            estimated_point_x_lf(point_id,ap_id) = vo(1,:);
            estimated_point_y_lf(point_id,ap_id) = vo(2,:);

        elseif (ap_id == 4) % 46

            estimated_point_x_lf(point_id,ap_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id)*(-1));
            estimated_point_y_lf(point_id,ap_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id));

            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x_lf(point_id,ap_id)));
            theta = deg2rad(-90);
            v = [estimated_point_x_lf(point_id,ap_id);estimated_point_y_lf(point_id,ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;
            vo = so + center; 
            estimated_point_x_lf(point_id,ap_id) = vo(1,:);
            estimated_point_y_lf(point_id,ap_id) = vo(2,:);
        end
        
         plot(estimated_point_x_lf(point_id, ap_id), estimated_point_y_lf(point_id, ap_id), "*")

    end
end


% Euclidean distance to ground truth, error
for ap_id=1:4
        
    
    errors = abs(sqrt( (estimated_point_x_lf(:, ap_id) - points_x').^2 + (estimated_point_y_lf(:, ap_id) - points_y').^2 ));

    figure;
    cdfplot(errors);
    title(['Location errors for AP ' num2str(ap_id) ' all points'])
end

save("mat_files/Coordinates/data_coordinates_lf_outdoor","estimated_point_x_lf", "estimated_point_y_lf");

cd(pwd_str)
