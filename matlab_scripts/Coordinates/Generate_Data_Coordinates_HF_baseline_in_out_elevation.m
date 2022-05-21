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
load("mat_files/in_out/HF/FTM/ftm_in_out.mat")
load('mat_files/in_out_map/data_distance_angle_true.mat');

estimated_distance = calculated_distance;

% load the angle
load("mat_files/in_out/HF/CSI/csi_in_out.mat")
angles_dp = -1*calculated_aoa;

% calibration
angles_dp(:,1,:) = angles_dp(:,1,:) - 1;
angles_dp(:,2,:) = angles_dp(:,2,:) - 2.2;


[n_points, routers, rotations] = size(calculated_aoa);

%% get the elevation
calculated_el = nan(size(calculated_aoa));
for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points
            el_aux = elevation_raw{id_point, id_router, id_rotation};
            if(~isempty(el_aux))
                calculated_el(id_point, id_router, id_rotation) = el_aux(1);
            end
        end
    end
end

offset = [+6,-2];
for ap = 1:routers
    calculated_el(:,ap,:) = calculated_el(:,ap,:) + offset(ap);
end

%% Localize: Estimate the positions

estimated_direct_path_angle    = angles_dp;% - median(abs(calculated_aoa - repmat(angles_RX, [1, 1, 8])), 'omitnan');
estimated_direct_path_distance = calculated_distance;

% https://en.wikipedia.org/wiki/Rotation_of_axes
for point_id=1:n_points
    for ap_id=1:routers
        for rotation_id=1:rotations
            
            estimated_point_x_p(point_id, ap_id, rotation_id) = sind(estimated_direct_path_angle(point_id,ap_id, rotation_id))*estimated_direct_path_distance(point_id,ap_id, rotation_id);
            estimated_point_y_p(point_id, ap_id, rotation_id) = cosd(estimated_direct_path_angle(point_id,ap_id, rotation_id))*estimated_direct_path_distance(point_id,ap_id, rotation_id);
            estimated_point_z_hf(id_point,id_router, id_rotation) = sind(calculated_el(id_point,id_router, id_rotation))*estimated_distance(id_point,id_router, id_rotation);

            % Rotate respect the AP
            x_center = RX_x(ap_id);
            y_center = RX_y(ap_id);

            if (ap_id == 1) % 43

                estimated_point_x_hf(point_id,ap_id, rotation_id) = RX_x(ap_id) + (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y_hf(point_id,ap_id, rotation_id) = RX_y(ap_id) + (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x_hf(point_id,ap_id, rotation_id)));
                theta = deg2rad(0);
                v = [estimated_point_x_hf(point_id,ap_id, rotation_id);estimated_point_y_hf(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x_hf(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y_hf(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 2) % 44

                estimated_point_x_hf(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y_hf(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x_hf(point_id,ap_id, rotation_id)));
                theta = deg2rad(-90);
                v = [estimated_point_x_hf(point_id,ap_id, rotation_id);estimated_point_y_hf(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;   
                vo = so + center; 
                estimated_point_x_hf(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y_hf(point_id,ap_id, rotation_id) = vo(2,:);

%             elseif (ap_id == 3) %45
% 
%                 estimated_point_x_hf(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
%                 estimated_point_y_hf(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));
% 
%                 % Now rotate the point
%                 center = repmat([x_center; y_center], 1, length(estimated_point_x_hf(point_id,ap_id, rotation_id)));
%                 theta = deg2rad(45);
%                 v = [estimated_point_x_hf(point_id,ap_id, rotation_id);estimated_point_y_hf(point_id,ap_id, rotation_id)];
%                 R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
%                 s = v - center;
%                 so = R*s;
%                 vo = so + center; 
%                 estimated_point_x_hf(point_id,ap_id, rotation_id) = vo(1,:);
%                 estimated_point_y_hf(point_id,ap_id, rotation_id) = vo(2,:);
% 
%             elseif (ap_id == 4) % 46
% 
%                 estimated_point_x_hf(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
%                 estimated_point_y_hf(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));
% 
%                 % Now rotate the point
%                 center = repmat([x_center; y_center], 1, length(estimated_point_x_hf(point_id,ap_id, rotation_id)));
%                 theta = deg2rad(-90);
%                 v = [estimated_point_x_hf(point_id,ap_id, rotation_id);estimated_point_y_hf(point_id,ap_id, rotation_id)];
%                 R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
%                 s = v - center;
%                 so = R*s;
%                 vo = so + center; 
%                 estimated_point_x_hf(point_id,ap_id, rotation_id) = vo(1,:);
%                 estimated_point_y_hf(point_id,ap_id, rotation_id) = vo(2,:);
            end
            
%             estimated_point_x_hf(point_id, ap_id, rotation_id);
%             estimated_point_y_hf(point_id, ap_id, rotation_id);
              plot(estimated_point_x_hf(point_id, ap_id, rotation_id), estimated_point_y_hf(point_id, ap_id, rotation_id), "*")

        end
    end
end

% Euclidean distance to ground truth, error
% errors = abs(sqrt( (estimated_point_x_hf(:, 2, 6) - points_x').^2 + (estimated_point_y_hf(:, 2, 6) - points_y').^2 ));
% figure;
% cdfplot(errors);

for ap_id=1:routers
        
    
    errors = abs(sqrt( (estimated_point_x_hf(:, ap_id, :) - points_x).^2 + (estimated_point_y_hf(:, ap_id, :) - points_y).^2 ));

    figure;
    cdfplot(errors);
    title(['Location errors for AP ' num2str(ap_id) ' all rotations'])
end

% for out metric
index_optimal_rotation = zeros(n_points,routers);

for ii = 1:n_points
    for jj = 1:routers
        [~, index_optimal_rotation(ii,jj)] = nanmin(estimated_distance(ii,jj,:));
    end
end
errors = zeros(n_points,routers);

for jj = 1:routers


    for ii = 1:n_points

        errors(ii,jj) = abs(sqrt( (estimated_point_x_hf(ii, jj, index_optimal_rotation(ii,jj)) - points_x(ii)).^2 + (estimated_point_y_hf(ii, jj, index_optimal_rotation(ii,jj)) - points_y(ii)).^2 ));

        
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

save("mat_files/Coordinates/data_coordinates_hf_in_out_elevation","estimated_point_x_hf", "estimated_point_y_hf", "estimated_point_z_hf");

cd(pwd_str)