clear 
close all
clc

pwd_str = pwd;
cd ../../


%% IMDEA testbed
load("mat_files/imdea/HF/CSI/csi_imdea.mat")
load("mat_files/imdea/HF/FTM/ftm_imdea.mat")
load("mat_files/imdea_map/data_distance_angle_true.mat")

[points, aps] = size(angles_RX);

index_optimal_ftm = zeros(size(angles_RX));
ftm_optimal = zeros(size(angles_RX));
aoa_optimal = zeros(size(angles_RX));

for ap = 1:3
    for point = 1:points
        
        calculated_distance_point = squeeze(calculated_distance(point,ap,:));
        [ftm_optimal(point,ap), index_min] = nanmin(calculated_distance_point);
        index_optimal_ftm(point,ap) = index_min;
        elevation_aux = elevation_raw{point,ap,index_min};
        if (~isempty(elevation_aux))
            aoa_optimal(point,ap) = elevation_raw{point,ap,index_min}(1);
        end

    end    
end

openfig("mat_files/imdea_map/Map.fig")

aoa_optimal(:,1) = aoa_optimal(:,1) - 2;
aoa_optimal(:,2) = aoa_optimal(:,2) - 6;
aoa_optimal(:,3) = aoa_optimal(:,3) + 2;

% aoa_errors = abs(aoa_optimal - angles_RX);
% ftm_errors = abs(ftm_optimal - distances_RX);
% 
figure, cdfplot((aoa_optimal(2:2:12,1)))
figure, cdfplot((aoa_optimal(:,2)))
figure, cdfplot((aoa_optimal(1:20,3)))
