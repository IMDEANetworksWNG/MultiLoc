clear all
clc
close all


%% General stuff
pwd_str = pwd;

cd ../../

addpath("functions")

load("mat_files/in_out_map/data_distance_angle_true.mat")
load("mat_files/in_out_map/points_coordinates.mat")
% load("mat_files/walls.mat")
load("mat_files/in_out_map/RX_coordinates.mat")
mkdir("mat_files/Coordinates/")


% RX_x = [RX_x_4, RX_x_5, RX_x_10];
% RX_y = [RX_y_4, RX_y_5, RX_y_10];



% angles_dp = angles_dp

% fig_map = Plot_Map_No_Points();


% load the ranging
load("mat_files/FTM/LF/estimated_distance_in_out.mat")
load('mat_files/in_out_map/data_distance_angle_true.mat');

% estimated_distance = calculated_distance;

% load the angle
load("mat_files/mD_track/LF/Data_3D_Direct_Path_in_out.mat")
% angles_dp = -1*angles_dp;


[n_points, routers] = size(angles_dp);

%% Localize: Estimate the positions

estimated_direct_path_angle    = angles_dp;% - median(abs(calculated_aoa - repmat(angles_RX, [1, 1, 8])), 'omitnan');
estimated_direct_path_distance = estimated_distance;

% https://en.wikipedia.org/wiki/Rotation_of_axes
for point_id=1:n_points
    for ap_id=1:routers
            
        estimated_point_x_p(point_id, ap_id) = sind(estimated_direct_path_angle(point_id,ap_id))*estimated_direct_path_distance(point_id,ap_id);
        estimated_point_y_p(point_id, ap_id) = cosd(estimated_direct_path_angle(point_id,ap_id))*estimated_direct_path_distance(point_id,ap_id);

        % Rotate respect the AP
        x_center = RX_x(ap_id);
        y_center = RX_y(ap_id);

        if (ap_id == 1) % 43

            estimated_point_x_lf(point_id,ap_id) = RX_x(ap_id) + (estimated_point_x_p(point_id,ap_id)*(-1));
            estimated_point_y_lf(point_id,ap_id) = RX_y(ap_id) + (estimated_point_y_p(point_id,ap_id));

            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x_lf(point_id,ap_id)));
            theta = deg2rad(0);
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
            theta = deg2rad(-90);
            v = [estimated_point_x_lf(point_id,ap_id);estimated_point_y_lf(point_id,ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;   
            vo = so + center; 
            estimated_point_x_lf(point_id,ap_id) = vo(1,:);
            estimated_point_y_lf(point_id,ap_id) = vo(2,:);

        end

    end
end


for jj = 1:routers

    errors = zeros(n_points,1);

    for ii = 1:n_points

        errors(ii,jj) = abs(sqrt( (estimated_point_x_lf(ii, jj) - points_x(ii)).^2 + (estimated_point_y_lf(ii, jj) - points_y(ii)).^2 ));

        
    end
    points = [];
    figure;
    if (jj == 1)
        points = 1:10;
    else
        points = 11:20;
    end
    cdfplot(errors(points,jj));
    title(['Location errors for AP ' num2str(jj) ' optimal rotation'])
end

save("mat_files/Coordinates/data_coordinates_lf_in_out","estimated_point_x_lf", "estimated_point_y_lf");

cd(pwd_str)