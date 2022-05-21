clear 
close all
clc

pwd_str = pwd;
cd ../../

load("mat_files/imdea/HF/CSI/angles_dp_amb.mat")
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
         aoa_optimal(point,ap) = angles_dp(point,ap, index_min);
    end    
end

openfig("mat_files/imdea_map/Map.fig")

aoa_optimal(:,1) = aoa_optimal(:,1) - 3;
aoa_optimal(:,3) = aoa_optimal(:,3) + 14;

aoa_errors = abs(aoa_optimal - angles_RX);
ftm_errors = abs(ftm_optimal - distances_RX);

figure, cdfplot(abs(aoa_errors(1:20,3)))
figure, cdfplot(abs(aoa_errors(2:2:12,1)))