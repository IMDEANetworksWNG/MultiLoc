clc
clear all
close all

%% Load the data
lf_aoa = load('mat_files/mD_track/LF/Data_3D_Direct_Path.mat');
hf_aoa = load('mat_files/indoor/HF/CSI/aoa_optimal_rotation.mat');
load('mat_files/indoor_map/Grid_LOS.mat')


hf_aoa = -1*hf_aoa.aoa_optimal_rotation;
lf_aoa = lf_aoa.angles_dp;

figure;
plot(hf_aoa(LOS_RX))
hold on
plot(lf_aoa(LOS_RX))
hold off

xlabel('Measuremnt')
ylabel('Degrees')
title('[LOS] Estimated angles for LF and HF')
legend('LF', 'HF')


difference = rmmissing(abs(hf_aoa(LOS_RX) - lf_aoa(LOS_RX)));

figure
cdfplot(difference)

%%
% What I need to do:
% For each point in LOS, check which of the rotations for HF best match the
% angle for LF
csi_hf = load('mat_files/indoor/HF/CSI/csi_indoor.mat');
ftm_hf = load('mat_files/indoor/HF/FTM/ftm_indoor.mat');
all_rotations_ftm = ftm_hf.calculated_distance;
all_rotations_aoa = csi_hf.calculated_aoa;
load('mat_files/indoor/HF/optimal_rotation_index.mat')
load('mat_files/indoor/HF/FTM/ftm_optimal_rotation.mat')
rotations = [0, 45, 270, 315, 90, 135, 180, 225];

chosen_rotation_based_on_lf = nan(110, 5);

% For each point
for ii=1:size(hf_aoa, 1)
    
    % For each AP
    for jj=1:size(hf_aoa, 2)

        lf_angle = lf_aoa(ii, jj);

        hf_angles = nan(8, 1);
        
        for kk=1:8
            
            hf_angles(kk) = -1*all_rotations_aoa(ii, jj, kk);
        end
        
        [minValue,closestIndex] = min(abs(hf_angles - lf_angle));
        
        if ~isnan(minValue)
            chosen_rotation_based_on_lf(ii, jj) = closestIndex;
        end
    end
end

% Generate a matrix with the estimated angles
estimated_angle_based_on_lf    = nan(110, 5);
estimated_distance_based_on_lf = nan(110, 5);

for ii=1:size(hf_aoa, 1)
    
    % For each AP
    for jj=1:size(hf_aoa, 2)
        
        estimated_rotation = chosen_rotation_based_on_lf(ii, jj);
        
        if ~isnan(estimated_rotation)
            estimated_angle_based_on_lf(ii, jj)    = all_rotations_aoa(ii, jj, estimated_rotation);
            estimated_distance_based_on_lf(ii, jj) = all_rotations_ftm(ii, jj, estimated_rotation);
        end
    end
end

% Load the ground truth
load('mat_files/indoor_map/data_distance_angle_true.mat')


figure;
%los
aux = abs(-1*hf_aoa(LOS_RX) - angles_RX(LOS_RX));
h = cdfplot(aux(:));
set(h, 'Color', 'r');
hold on
aux = abs(estimated_angle_based_on_lf(LOS_RX) - angles_RX(LOS_RX));
h = cdfplot(aux(:));
set(h, 'Color', 'b');
% nlos
aux = abs(-1*hf_aoa(~LOS_RX) - angles_RX(~LOS_RX));
h = cdfplot(aux(:));
set(h,'Marker', '+', 'Color', 'r');
hold on
aux = abs(estimated_angle_based_on_lf(~LOS_RX) - angles_RX(~LOS_RX));
h = cdfplot(aux(:));
set(h,'Marker', '+', 'Color', 'b');
hold off

xlabel('Error in degrees')
title('[Indoor AoA] Optimal rotation vs rotation taken using LF')
legend('[LOS] Optimal rotation', '[LOS] Rotation asisted by LF', '[NLOS] Optimal rotation', '[NLOS] Rotation asisted by LF')



% Same for ToF
figure;
% los
aux = abs(ftm_optimal_rotation(LOS_RX) - distances_RX(LOS_RX));
h = cdfplot(aux(:));
set(h, 'Color', 'r');
hold on
aux = abs(estimated_distance_based_on_lf(LOS_RX) - distances_RX(LOS_RX));
h = cdfplot(aux(:));
set(h, 'Color', 'b');
% nlos
aux = abs(ftm_optimal_rotation(~LOS_RX) - distances_RX(~LOS_RX));
h = cdfplot(aux(:));
set(h,'Marker', '+', 'Color', 'r');
hold on
aux = abs(estimated_distance_based_on_lf(~LOS_RX) - distances_RX(~LOS_RX));
h = cdfplot(aux(:));
set(h,'Marker', '+', 'Color', 'b');
hold off

xlabel('Error in meters')
title('[Indoor ToF] Optimal rotation vs rotation taken using LF')
legend('[LOS] Optimal rotation', '[LOS] Rotation asisted by LF', '[NLOS] Optimal rotation', '[NLOS] Rotation asisted by LF')
