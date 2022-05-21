clear all
clc
close all


%% General stuff
pwd_str = pwd;

cd ../../

addpath("functions")

load("mat_files/indoor_map/data_distance_angle_true.mat")
load("mat_files/indoor_map/points_coordinates.mat")
% load("mat_files/walls.mat")
load("mat_files/indoor_map/RX_coordinates.mat")
mkdir("mat_files/LOS_NLOS/")

load('mat_files/indoor_map/Grid_LOS.mat')


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

% load the angle
load("mat_files/indoor/HF/CSI/csi_indoor.mat")
angles_dp = calculated_aoa;
% load the offset
load("mat_files/indoor/HF/CSI/aoa_offset.mat")

std_all = zeros(size(azimut_raw));
los_factor_all = nan(size(azimut_raw));

[n_points, routers, rotations] = size(calculated_aoa);


for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points
            mag = magnitudes_raw(id_point,id_router ,id_rotation);
            mag = mag{1,1};
            std_all(id_point, id_router, id_rotation) = median(std(mag./mean(mag)));
%             std_all(id_point, id_router, id_rotation) = median(median(mag));
            power_aux = power_raw{id_point, id_router, id_rotation};
            if (~isempty(power_aux))
                los_factor_all(id_point, id_router, id_rotation) = max(power_aux)/sum(power_aux);
            end
        end
    end
end

LOS_RX = repmat(LOS_RX,1,1,8);
los_factor_all_los = los_factor_all(LOS_RX);
los_factor_all_nlos = los_factor_all(~LOS_RX);

% LOS_RX = repmat(LOS_RX,1,1,8);
std_los = std_all(LOS_RX);
std_nlos = std_all(~LOS_RX);

figure, cdfplot(std_los)
hold on
cdfplot(std_nlos)
legend("LOS", "NLOS")

figure, cdfplot(los_factor_all_los)
hold on
cdfplot(los_factor_all_nlos)
legend("LOS", "NLOS")

std_los = std_all.* LOS_RX;
std_los_nans = std_los;
std_los_nans(~LOS_RX) = nan;

std_nlos = std_all.* (~LOS_RX);
std_nlos_nans = std_nlos;
std_nlos_nans(LOS_RX) = nan;

save("mat_files/LOS_NLOS/std_all", "std_all");

% for id_rotation = 1:8
%     figure, cdfplot(reshape(squeeze(std_los_nans(:,:,id_rotation)),n_points*routers,1))
%     hold on
%     cdfplot(reshape(squeeze(std_nlos_nans(:,:,id_rotation)),n_points*routers,1))
%     legend("LOS", "NLOS") 
%     
% end
% mag_los = magnitudes_raw(20,5,3);
% mag_los = mag_los{1,1};

% mag_nlos = magnitudes_raw(27,5,1);
% mag_nlos = mag_nlos{1,1};

% std_los = median(std(mag_los./mean(mag_los)))
% std_nlos = median(std(mag_nlos./mean(mag_nlos)))

cd(pwd_str)