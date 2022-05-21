clear all
close all
clc

%% Indoor scenario
csi_indoor  = load('mat_files/indoor/HF_15db_mD-Track/CSI/csi_indoor.mat');
ftm_indoor = load('mat_files/indoor/HF_15db_mD-Track/FTM/ftm_indoor.mat');

csi_indoor_15  = load('mat_files/indoor/HF_15db_mD-Track/CSI/csi_indoor.mat');
ftm_indoor_15 = load('mat_files/indoor/HF_15db_mD-Track/FTM/ftm_indoor.mat');

% Map data
indoor_rx                   = load('Maps/Map_indoor/RX_coordinates.mat');
indoor_coordinates          = load('Maps/Map_indoor/points_coordinates.mat');
indoor_true_angles_distance = load('Maps/Map_indoor/data_distance_angle_true.mat');

load('mat_files/indoor/HF/optimal_rotation_index.mat')

% Indexes for LOS/NLOS
load('mat_files/indoor_map/Grid_LOS.mat')

los_grid = repmat(LOS_RX, [1 1 8]);

% LOS factor
los_factor_los       = nan(110, 5, 8);
number_of_paths_los  = nan(110, 5, 8);
los_factor_nlos      = nan(110, 5, 8);
number_of_paths_nlos = nan(110, 5, 8);

for ii=1:110
    
    % For each AP
    for jj=1:5

        % For each rotation
        for kk=1:8
            
            power = csi_indoor.power_raw{ii, jj, kk};
            
            if size(power, 2) > 0
                
                if LOS_RX(ii, jj)
                    los_factor_los(ii, jj, kk)      = max(power)/sum(power);
                    number_of_paths_los(ii, jj, kk) = size(power, 2);                
                else
                    los_factor_nlos(ii, jj, kk)      = max(power)/sum(power);
                    number_of_paths_nlos(ii, jj, kk) = size(power, 2);         
                end
            end
        end
    end
end

% figure;
% los = los_factor(LOS_RX);
% nlos = los_factor(~LOS_RX);
% group = [    ones(size(los));
%          2 * ones(size(nlos));
%          ];
% 
% boxplot([los; nlos], group);
% set(gca,'XTickLabel',{'LOS','NLOS'})


los = los_factor_los;
nlos = los_factor_nlos;

% Los factor per rotation
figure;
hold on
for ii=1:8
    aux = los(:, :, ii);
    cdfplot(aux(:));
end
hold off
legend(string(0:45:315))
xlabel('LOS Factor per rotation')
title('[LOS] LOS Factor')

figure;
hold on
for ii=1:8
    aux = nlos(:, :, ii);
    cdfplot(aux(:));
end
hold off
legend(string(0:45:315))
xlabel('LOS Factor per rotation')
title('[NLOS] LOS Factor')

figure;
cdfplot(los(:));
hold on
cdfplot(nlos(:));

%set(h,'Marker', 'o');

legend('LOS', 'NLOS')
xlabel('LOS Factor')
title('LOS Factor')

% Number of paths
los = number_of_paths_los;
nlos = number_of_paths_nlos;

figure;
cdfplot(los(:));
hold on
cdfplot(nlos(:));
legend('LOS', 'NLOS')
xlabel('Number of paths')
title('Number of paths')
