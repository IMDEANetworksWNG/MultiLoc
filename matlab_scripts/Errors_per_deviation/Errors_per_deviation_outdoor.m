clc
clear all
close all


cd ../../
%% Indoor scenario
csi_outdoor  = load('mat_files/outdoor/HF/CSI/csi_outdoor.mat');
ftm_outdoor  = load('mat_files/outdoor/HF/FTM/ftm_outdoor.mat');

% Map data
outdoor_rx                  = load('mat_files/outdoor_map/RX_coordinates.mat');
outoor_coordinates          = load('mat_files/outdoor_map/points_coordinates.mat');
outoor_true_angles_distance = load('mat_files/outdoor_map/data_distance_angle_true.mat');

rgb = viridis(9);

load('mat_files/outdoor/HF/optimal_rotation_index.mat')

% Calibration
calibration = load('mat_files/outdoor/HF/CSI/calibration.mat');

aoa_optimal_rotation = nan(16, 4);
ftm_optimal_rotation = nan(16, 4);

% Save the distance and FTM for the optimal rotation
for ii=1:size(outoor_coordinates.points_x, 2)
    
    % For each AP
    for jj=1:size(outdoor_rx.RX_x, 2)

        kk = optimal_rotation_index(ii, jj);
        aoa_optimal_rotation(ii, jj) = csi_outdoor.calculated_aoa(ii, jj, kk);
        ftm_optimal_rotation(ii, jj) = ftm_outdoor.calculated_distance(ii, jj, kk);
        
    end
end

% We need to have the true_angles and true_distances copied for each
% rotation
true_angles    = repmat(outoor_true_angles_distance.angles_RX, [1 1 8]);
true_distances = repmat(outoor_true_angles_distance.distances_RX, [1 1 8]);

rotations = sort([0, 45, 270, 315, 90, 135, 180, 225]);

errors_per_deviation = nan(size(outoor_coordinates.points_x, 2), size(outdoor_rx.RX_x, 2), size(rotations, 2));

% For each point
for ii=1:size(outoor_coordinates.points_x, 2)
    
    % For each AP
    for jj=1:size(outdoor_rx.RX_x, 2)

        % Angle for that point and AP for the optimal rotation
        optimal_angle_index = optimal_rotation_index(ii, jj);
        optimal_angle_gt = rotations(optimal_angle_index);
        
        for rotation_index=1:size(rotations, 2)
            
            current_rotation = rotations(rotation_index);
            offset           = mod(optimal_angle_gt - current_rotation, 360);
            
            save_index = rotations==offset;
            
            errors_per_deviation(ii, jj, save_index) = abs(csi_outdoor.calculated_aoa(ii, jj, rotations==current_rotation) - -1*true_angles(ii, jj, rotations==current_rotation));
        end
    end
end

% Calibrate
errors_per_deviation = abs(errors_per_deviation - calibration.error_aoa);

% Remove outliers
errors_per_deviation(errors_per_deviation>150) = nan;

figure
hold on

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

for rotation_index=1:size(rotations, 2)
    
    aux = errors_per_deviation(:, :, rotation_index);
    h = cdfplot(aux(:));
    
    set(h,'Marker', markers{rotation_index}, 'LineStyle', '-', 'color', rgb(rotation_index, :), 'LineWidth', 1);
end

% LF AoA
load('mat_files/mD_track/LF/Data_3D_Direct_Path_outdoor.mat')
angles_dp(:, [1 4]) = angles_dp(:, [4 1]);

% Calibrate
angles_dp(:,2) = angles_dp(:,2) + 3.5;
angles_dp(:,3) = angles_dp(:,3) - 2;

error = abs(outoor_true_angles_distance.angles_RX - angles_dp);
h = cdfplot(error(:));
set(h, 'LineStyle', '-', 'color', rgb(9, :), 'LineWidth', 2);

hold off
xlabel('Error in degrees')
title('[Outdoor] AoA errors if we chose the wrong rotation by N degrees')
legend([string(rotations) 'LF'], 'Location', "SouthEast")
matlab2tikz('plots/final_plots/rotation_comparison/AoA_LOS_Outdoor.tikz');


% FTM for LOS and NLOS
errors_per_deviation  = nan(size(outoor_coordinates.points_x, 2), size(outdoor_rx.RX_x, 2), size(rotations, 2));

% For each point
for ii=1:size(outoor_coordinates.points_x, 2)
    
    % For each AP
    for jj=1:size(outdoor_rx.RX_x, 2)

        % Angle for that point and AP for the optimal rotation
        optimal_angle_index = optimal_rotation_index(ii, jj);
        optimal_angle_gt = rotations(optimal_angle_index);
        
        for rotation_index=1:size(rotations, 2)
            
            current_rotation = rotations(rotation_index);
            offset           = mod(optimal_angle_gt - current_rotation, 360);
            
            save_index = rotations==offset;
            
            errors_per_deviation(ii, jj, save_index)  = abs(ftm_outdoor.calculated_distance(ii, jj, rotations==current_rotation) - true_distances(ii, jj, rotations==current_rotation));
        end
    end
end

figure
hold on

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

for rotation_index=1:size(rotations, 2)
    
    aux = errors_per_deviation(:, :, rotation_index);
    h = cdfplot(aux(:));
    
    set(h,'Marker', markers{rotation_index}, 'LineStyle', '-', 'color', rgb(rotation_index, :), 'LineWidth', 1);
end

% LF FTM outdoor
load('mat_files/FTM/LF/estimated_distance_outdoor.mat')

estimated_distance(:, [1 4]) = estimated_distance(:, [4 1]);

% LF FTM
error = abs(outoor_true_angles_distance.distances_RX - estimated_distance);
h = cdfplot(error(:));
set(h, 'LineStyle', '-', 'color', rgb(9, :), 'LineWidth', 2);

hold off
xlabel('Error in meters')
title('[Outdoor LOS] ToF errors if we chose the wrong rotation by N degrees')
legend([string(rotations) 'LF'], 'Location', "SouthEast")

matlab2tikz('plots/final_plots/rotation_comparison/ToF_LOS_Outdoor.tikz');




%% Outdoor scenario
% csi_outdoor = load('mat_files/outdoor/HF/CSI/csi_outdoor.mat');
% csi_indoor  = load('mat_files/outdoor/HF/FTM/ftm_outdoor.mat');
% 
% % Map data
% outdoor_rx                    = load('mat_files/outdoor_map/RX_coordinates.mat');
% outdoor_coordinates           = load('mat_files/outdoor_map/points_coordinates.mat');
% outdoor_true_angles_distances = load('mat_files/outdoor_map/data_distance_angle_true.mat');
% 
% % We need to have the true_angles and true_distances copied for each
% % rotation
% true_angles    = repmat(outdoor_true_angles_distances.angles_RX, [1 1 8]);
% true_distances = repmat(outdoor_true_angles_distances.distances_RX, [1 1 8]);
% 
% % Angle for each point, rotation and AP
% figure
% cdfplot(abs(csi_outdoor.calculated_aoa(:) - true_angles(:)));
% xlabel('Error in degrees')
% title('[Outdoor] AoA errors for all points, all APs, all rotations')