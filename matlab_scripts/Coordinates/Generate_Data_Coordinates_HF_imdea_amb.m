clear all
clc
close all


%% General stuff
pwd_str = pwd;

cd ../../

addpath("functions")

load("mat_files/imdea_map/data_distance_angle_true.mat")
load("mat_files/imdea_map/points_coordinates.mat")
% load("mat_files/walls.mat")
load("mat_files/imdea_map/RX_coordinates.mat")
mkdir("mat_files/Coordinates/")


% RX_x = [RX_x_4, RX_x_5, RX_x_10];
% RX_y = [RX_y_4, RX_y_5, RX_y_10];



% angles_dp = angles_dp

% fig_map = Plot_Map_No_Points();


% load the ranging
load("mat_files/imdea/HF/FTM/ftm_imdea.mat")
load('mat_files/imdea_map/data_distance_angle_true.mat');

estimated_distance = calculated_distance;

% load the angle
load("mat_files/mD_track/LF/Data_3D_Direct_Path_imdea.mat")
angles_dp_lf = angles_dp;
angles_dp_lf(:,3) = angles_dp_lf(:,3) - 6;
angles_dp_lf(:,2) = angles_dp_lf(:,2) + 6;
angles_dp_lf([19,20,27:31],1) = angles_dp_lf([19,20,27:31],1)*(-1);

% load the angle
load("mat_files/imdea/HF/CSI/csi_imdea.mat")
angles_dp = -1*calculated_aoa;

% calibration
cal = [3, 0, 14];
% angles_dp(:,1,:) = angles_dp(:,1,:) + 3;
% angles_dp(:,3,:) = angles_dp(:,3,:) - 14;


[n_points, routers, rotations] = size(calculated_aoa);

%% remove ambiguity
angles_dp = calculated_aoa;
angles_dp = angles_dp*(-1);
% angles_dp_2 = angles_dp;

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
    
    
    [row, column, ~] = ind2sub(size(angles_dp), ii);
    
    possible_angles = [angles_dp_aux angles_dp_aux_amb] + cal(column);

    
    angle_dp_lf_aux = angles_dp_lf(row, column);
    
    [~, index_min] = min(abs(angle_dp_lf_aux - possible_angles));
    angles_dp(ii) = possible_angles(index_min);
end
angles_dp = angles_dp;
% save("angles_dp"
save("mat_files/imdea/HF/CSI/angles_dp_amb.mat","angles_dp")

angles_dp(:,1,:) = angles_dp(:,1,:) + 3;
angles_dp(:,3,:) = angles_dp(:,3,:) - 14;

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

offset = [-2, -8, +2];
for ap = 1:routers
    calculated_el(:,ap,:) = calculated_el(:,ap,:) + offset(ap);
end

angles_dp = real(asind(sind(angles_dp)./cosd(calculated_el)));

%% Localize: Estimate the positions

estimated_direct_path_angle    = angles_dp;% - median(abs(calculated_aoa - repmat(angles_RX, [1, 1, 8])), 'omitnan');
estimated_direct_path_distance = calculated_distance;

% https://en.wikipedia.org/wiki/Rotation_of_axes
for point_id=1:n_points
    for ap_id=1:routers
        for rotation_id=1:rotations
            
            estimated_point_x_p(point_id, ap_id, rotation_id) = sind(estimated_direct_path_angle(point_id,ap_id, rotation_id))*estimated_direct_path_distance(point_id,ap_id, rotation_id);
            estimated_point_y_p(point_id, ap_id, rotation_id) = cosd(estimated_direct_path_angle(point_id,ap_id, rotation_id))*estimated_direct_path_distance(point_id,ap_id, rotation_id);

            % Rotate respect the AP
            x_center = RX_x(ap_id);
            y_center = RX_y(ap_id);

            if (ap_id == 1) % 43

                estimated_point_x_hf(point_id,ap_id, rotation_id) = RX_x(ap_id) + (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y_hf(point_id,ap_id, rotation_id) = RX_y(ap_id) + (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x_hf(point_id,ap_id, rotation_id)));
                theta = deg2rad(180);
                v = [estimated_point_x_hf(point_id,ap_id, rotation_id);estimated_point_y_hf(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x_hf(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y_hf(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 2) % 52

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

            elseif (ap_id == 3) % 54

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
    if (jj == 3)
        points = 1:10;
    elseif (jj == 2)
        points = 21:26;
    else
        points = 2:2:12;
    end
    cdfplot(errors(points,jj));
    title(['Location errors for AP ' num2str(jj) ' optimal rotation'])
end

save("mat_files/Coordinates/data_coordinates_hf_imdea","estimated_point_x_hf", "estimated_point_y_hf");

cd(pwd_str)