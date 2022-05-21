clear all
clc
close all


%% General stuff
%pwd_str = pwd;

%cd ../../

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
estimated_distance = calculated_distance;

% load the angle
load("mat_files/outdoor/HF/CSI/csi_outdoor.mat")
angles_dp = calculated_aoa;

routers = 4;
n_points = 16;

[n_points, routers, rotations] = size(calculated_aoa);

%% Estimate the positions
% estimated_point_x_hf = zeros(size(calculated_aoa));
% estimated_point_y_hf = zeros(size(calculated_aoa));
% 
% for id_rotation = 1:rotations
%     for id_router = 1:routers
%         for id_point = 1:n_points
% 
%             estimated_point_x_aux = sind(angles_dp(id_point,id_router, id_rotation))*estimated_distance(id_point,id_router, id_rotation);
%             estimated_point_y_aux = cosd(angles_dp(id_point,id_router, id_rotation))*estimated_distance(id_point,id_router, id_rotation);
% 
% 
%     %         figure, plot(estimated_point_x_p(id_point,id_router), estimated_point_y_p(id_point,id_router), "*")
% 
% 
%             if (id_router == 1 || id_router == 3 || id_router == 4)
%                 estimated_point_x_hf(id_point,id_router, id_rotation) = RX_x(id_router) - (estimated_point_x_aux);
%                 estimated_point_y_hf(id_point,id_router, id_rotation) = RX_y(id_router) - (estimated_point_y_aux*(-1));
%             else
%                 estimated_point_x_hf(id_point,id_router, id_rotation) = RX_x(id_router) - (estimated_point_x_aux*(-1));
%                 estimated_point_y_hf(id_point,id_router, id_rotation) = RX_y(id_router) - estimated_point_y_aux;
%             end
%         end
%     end
% end
% 
% % Euclidean distance to ground truth, error
% errors = abs(sqrt( (estimated_point_x_hf(:, 2, 6) - points_x').^2 + (estimated_point_y_hf(:, 2, 6) - points_y').^2 ));
% figure;
% cdfplot(errors);

%% Localize: Estimate the positions

estimated_direct_path_angle    = calculated_aoa - median(angdiff(calculated_aoa, repmat(angles_RX, [1, 1, 8])), 'omitnan');
estimated_direct_path_distance = calculated_distance;

% AoA ToF errors
% for ap_id=1
%     
%     aoa_errors = abs(angdiff(estimated_direct_path_angle(:, ap_id), angles_RX(:, ap_id)));
%     tof_errors = abs(estimated_direct_path_distance(:, ap_id) - distances_RX(:, ap_id));
% 
%     figure
%     subplot(1,2,1);
%     cdfplot(aoa_errors)
%     title(['AoA error ap ' num2str(ap_id)])
% 
%     subplot(1,2,2);
%     cdfplot(tof_errors)
%     title(['ToF error ' num2str(ap_id)])
% 
%     eval(['aoa_errors_'  num2str(ap_id) ' = aoa_errors;']);
%     eval(['tof_errors_'  num2str(ap_id) ' = tof_errors;']);
%     
% end

% Remove the median

% uiopen('/home/imdea/Documents/MATLAB/Music_Mikrotik/mat_files/outdoor_map/outdoor_map.fig',1)
% hold on
% pause(1)

% https://en.wikipedia.org/wiki/Rotation_of_axes
for point_id=1:16
    for ap_id=1:4
        for rotation_id=1:8
            
            estimated_point_x_p(point_id, ap_id, rotation_id) = sind(estimated_direct_path_angle(point_id,ap_id, rotation_id)*-1)*estimated_direct_path_distance(point_id,ap_id, rotation_id);
            estimated_point_y_p(point_id, ap_id, rotation_id) = cosd(estimated_direct_path_angle(point_id,ap_id, rotation_id)*-1)*estimated_direct_path_distance(point_id,ap_id, rotation_id);

            % Rotate respect the AP
            x_center = RX_x(ap_id);
            y_center = RX_y(ap_id);

            if (ap_id == 1) % 43

                estimated_point_x(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id));
                estimated_point_y(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id)*(-1));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_id, rotation_id)));
                theta = deg2rad(45);
                v = [estimated_point_x(point_id,ap_id, rotation_id);estimated_point_y(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 2) % 44

                estimated_point_x(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_id, rotation_id)));
                theta = deg2rad(-45);
                v = [estimated_point_x(point_id,ap_id, rotation_id);estimated_point_y(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;   
                vo = so + center; 
                estimated_point_x(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 3) %45

                estimated_point_x(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_id, rotation_id)));
                theta = deg2rad(35);
                v = [estimated_point_x(point_id,ap_id, rotation_id);estimated_point_y(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y(point_id,ap_id, rotation_id) = vo(2,:);

            elseif (ap_id == 4) % 46

                estimated_point_x(point_id,ap_id, rotation_id) = RX_x(ap_id) - (estimated_point_x_p(point_id,ap_id, rotation_id)*(-1));
                estimated_point_y(point_id,ap_id, rotation_id) = RX_y(ap_id) - (estimated_point_y_p(point_id,ap_id, rotation_id));

                % Now rotate the point
                center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_id, rotation_id)));
                theta = deg2rad(-90);
                v = [estimated_point_x(point_id,ap_id, rotation_id);estimated_point_y(point_id,ap_id, rotation_id)];
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                s = v - center;
                so = R*s;
                vo = so + center; 
                estimated_point_x(point_id,ap_id, rotation_id) = vo(1,:);
                estimated_point_y(point_id,ap_id, rotation_id) = vo(2,:);
            end
            
%             estimated_point_x(point_id, ap_id, rotation_id);
%             estimated_point_y(point_id, ap_id, rotation_id);
%              plot(estimated_point_x(point_id, ap_id, rotation_id), estimated_point_y(point_id, ap_id, rotation_id), "*")

        end
    end
end

% Euclidean distance to ground truth, error
% errors = abs(sqrt( (estimated_point_x(:, 2, 6) - points_x').^2 + (estimated_point_y(:, 2, 6) - points_y').^2 ));
% figure;
% cdfplot(errors);

for ap_id=1:4
        
    
    errors = abs(sqrt( (estimated_point_x(:, ap_id, :) - points_x').^2 + (estimated_point_y(:, ap_id, :) - points_y').^2 ));

    figure;
    cdfplot(errors);
    title(['Location errors for AP ' num2str(ap_id) ' all rotations'])
end

% save("mat_files/Coordinates/data_coordinates_hf_outdoor","estimated_point_x", "estimated_point_y");

%cd(pwd_str)