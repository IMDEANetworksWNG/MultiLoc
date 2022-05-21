clear 
close all
clc

%%
pwd_str = pwd;

cd ../../../

mkdir("plots")
mkdir("plots/final_plots_subplot_new")

addpath('functions/');
% addpath('matlab_scripts/fu');

scenarios = ["indoor", "outdoor", "indoor", "in_out", "imdea"];
scenarios_save = ["indoor", "outdoor", "2ap", "in_out", "imdea"];

subarrays = [8,4,2,1];
index_subarrays = {0:7,0:2:6,0:4:4,0};


colors = 	[0, 0.4470, 0.7410,
          	0.8500, 0.3250, 0.0980,
          	0.9290, 0.6940, 0.1250,
          	0.4940, 0.1840, 0.5560,
          	0.4660, 0.6740, 0.1880];

mkdir("mat_files/")
mkdir("mat_files/th_angle_distance")
        
th_disagreement_all = [4.5,4.5,2,2];


% for scenario = 1:length(scenarios)
for scenario = 1
% for scenario = 2
mkdir(strcat("mat_files/CDF_data_new/", scenarios_save(scenario)))
% load the coordinates of the client
load(strcat("mat_files/",scenarios(scenario),"_map/points_coordinates.mat"))
load(strcat("mat_files/",scenarios(scenario),"_map/data_distance_angle_true.mat"))

% load the HF FTM data
load(strcat("mat_files/",scenarios(scenario),"/HF/FTM/ftm_",scenarios(scenario),".mat"))
load(strcat("mat_files/",scenarios(scenario),"/HF/CSI/csi_",scenarios(scenario),".mat"))

if (scenario == 1)
    load(strcat("mat_files/",scenarios(scenario),"/HF/calculated_aoa_amb.mat"))
    calculated_aoa = angles_dp;
end
% calculated_aoa(:,1,:) = calculated_aoa(:,1,:) - aoa_offset(1);
% calculated_aoa(:,2,:) = calculated_aoa(:,2,:) - aoa_offset(2);
% calculated_aoa(:,3,:) = calculated_aoa(:,3,:) - aoa_offset(3);
% calculated_aoa(:,4,:) = calculated_aoa(:,4,:) - aoa_offset(4);
% calculated_aoa(:,5,:) = calculated_aoa(:,5,:) - aoa_offset(5);


% load the LF FTM data
load(strcat("mat_files/FTM/LF/estimated_distance_",scenarios(scenario),".mat"))

load(strcat("mat_files/mD_track/LF/Data_3D_Direct_Path_",scenarios(scenario),".mat"))
angles_dp(1:74,2) = angles_dp(1:74,2) - 3.5;
angles_dp(26:38,5) = angles_dp(26:38,5) - 5;
angles_dp(1:74,3) = angles_dp(1:74,3) - 7;



% load the HF data
load(strcat("mat_files/Coordinates/data_coordinates_hf_",scenarios(scenario),".mat"))
% load the LF data
load(strcat("mat_files/Coordinates/data_coordinates_lf_",scenarios(scenario),".mat"))

% load the multiband baseline
load(strcat("mat_files/Coordinates/data_coordinates_multiband_baseline_",scenarios(scenario),".mat"))

% load the hf baseline
load(strcat("mat_files/Coordinates/data_coordinates_hf_bs_",scenarios(scenario),".mat"))


[n_points, routers, rotations] = size(estimated_point_x_hf);

% fig = figure;

if (scenario == 2)
    points_x = points_x.';
    points_y = points_y.';
    estimated_distance(:, [1 4]) = estimated_distance(:, [4 1]);

elseif (scenario == 3)
    router_indexes = [2 5];
    calculated_distance = calculated_distance(:, router_indexes, :);
    estimated_distance = estimated_distance(:, router_indexes);
    angles_dp = angles_dp(:, router_indexes);

    estimated_point_x_hf = estimated_point_x_hf(:, router_indexes, :);
    estimated_point_y_hf = estimated_point_y_hf(:, router_indexes, :);
    estimated_point_x_lf = estimated_point_x_lf(:, router_indexes, :);
    estimated_point_y_lf = estimated_point_y_lf(:, router_indexes, :);

    estimated_point_x_mb_bs = estimated_point_x_mb_bs(:, router_indexes, :);
    estimated_point_y_mb_bs = estimated_point_y_mb_bs(:, router_indexes, :);

    estimated_point_x_hf_bs = estimated_point_x_hf_bs(:, router_indexes, :);
    estimated_point_y_hf_bs = estimated_point_y_hf_bs(:, router_indexes, :);
    [n_points, routers, rotations] = size(estimated_point_x_hf);
    
end
alpha = 0.85;
% for ii = 1:length(subarrays)
    
for ii = 1:length(subarrays)   
estimated_point_x_hf_rotations = nan(n_points,routers,rotations/subarrays(ii));
estimated_point_y_hf_rotations = nan(n_points,routers,rotations/subarrays(ii));

disagreement = zeros(n_points,routers,rotations/subarrays(ii));
chosen_rotation_with_hf = nan(n_points,routers, rotations/subarrays(ii));

optimal_aoa = nan(size(estimated_point_x_hf_rotations));
optimal_ftm = nan(size(estimated_point_x_hf_rotations));


for id_point = 1:n_points
    for id_router = 1:routers
        for id_rotation = 1:(rotations/subarrays(ii))
            
            % take ii
            estimated_point_x_hf_aux = estimated_point_x_hf(id_point, id_router, id_rotation + index_subarrays{ii});
            estimated_point_y_hf_aux = estimated_point_y_hf(id_point, id_router, id_rotation + index_subarrays{ii});
            
            calculated_aoa_aux = calculated_aoa(id_point, id_router, id_rotation + index_subarrays{ii});
            calculated_distance_aux = calculated_distance(id_point, id_router, id_rotation + index_subarrays{ii});
            
            % take the lf one
            estimated_point_x_lf_aux = estimated_point_x_lf(id_point, id_router);
            estimated_point_y_lf_aux = estimated_point_y_lf(id_point, id_router);
            
            % take the one that better agree with LF
            distance_arrays = sqrt((estimated_point_x_hf_aux - estimated_point_x_lf_aux).^2 + (estimated_point_y_hf_aux - estimated_point_y_lf_aux).^2);
            
            % check if we have coverage, if not take the lf information
            if (sum(isnan(distance_arrays)) == subarrays(ii))
                disagreement(id_point, id_router, id_rotation) = Inf;
            else
                % if we have coverage, take the subarray that better agrees
                % with LF
%                 [~, index_min_rotation] = min(distance_arrays);
                [~, index_min_rotation] = min(squeeze(calculated_distance(id_point, id_router,id_rotation + index_subarrays{ii})));
                estimated_point_x_hf_aux = estimated_point_x_hf_aux(index_min_rotation);
                estimated_point_y_hf_aux = estimated_point_y_hf_aux(index_min_rotation);
                distance_lf_hf = sqrt((estimated_point_x_hf_aux - estimated_point_x_lf_aux)^2 + (estimated_point_y_hf_aux - estimated_point_y_lf_aux)^2);
                disagreement(id_point, id_router, id_rotation) = distance_lf_hf;
                chosen_rotation_with_hf(id_point, id_router, id_rotation) = index_min_rotation;
                estimated_point_x_hf_rotations(id_point, id_router, id_rotation) = estimated_point_x_hf_aux;
                estimated_point_y_hf_rotations(id_point, id_router, id_rotation) = estimated_point_y_hf_aux;
                optimal_aoa(id_point, id_router, id_rotation) = calculated_aoa_aux(index_min_rotation);
                optimal_ftm(id_point, id_router, id_rotation) = calculated_distance_aux(index_min_rotation);

            end
            
        end
    end
end

%% Get the ap chosen by multiloc
index_ap_in = false( n_points, routers, rotations/subarrays(ii));

for id_rotation = 1:(rotations/subarrays(ii))
        
%     estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
%     estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
%     estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    index_nan = isnan(squeeze(estimated_point_y_hf_rotations(:,:,id_rotation)));
%     estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
    distances = estimated_distance;
    distances(index_nan) = nan;
    
    min_distance = nanmin(distances,[],2);
    index_in_disagreement = min_distance == distances;
    index_ap_in(:,:,id_rotation) = index_in_disagreement;
end

%% Get the errors

% figure,
error_all = [];

% compute the errors
error_hf = sqrt((estimated_point_x_hf_rotations - points_x).^2 + (estimated_point_y_hf_rotations - points_y).^2);
error_hf(~index_ap_in) = nan;
error_hf = error_hf(:);

estimated_point_x_lf_rot = repmat(estimated_point_x_lf,1,1,rotations/subarrays(ii));
estimated_point_y_lf_rot = repmat(estimated_point_y_lf,1,1,rotations/subarrays(ii));

error_lf = sqrt((estimated_point_x_lf_rot - points_x).^2 + (estimated_point_y_lf_rot - points_y).^2);
error_lf(~index_ap_in) = nan;
error_lf = error_lf(:);

index_nan = isnan(error_hf);
error_lf(index_nan) = nan;

error_all(:,1) = error_hf;
error_all(:,2) = error_lf;

[~, index_min] = nanmin(error_all,[],2);
index_min(index_nan) = nan;

in_hf = index_min == 1;
in_lf = index_min == 2;

%% angle
angles_dp_aux = repmat(angles_dp,1,1,rotations/subarrays(ii));
% optimal_aoa = optimal_aoa*(-1);
error_angle = angles_dp_aux - optimal_aoa;
error_angle = abs(error_angle(:));
error_angle(~index_ap_in) = nan;

error_angle_hf = error_angle;
error_angle_hf(~in_hf) = nan;

error_angle_lf = error_angle;
error_angle_lf(~in_lf) = nan;


figure,
subplot(3,1,1)
histogram(error_angle_hf(:), 100,'Normalization','probability')
hold on
histogram(error_angle_lf(:), 100,'Normalization','probability')
title("Angle difference")
legend("Better HF", "Better LF")

%% Distance
estimated_distance_aux = repmat(estimated_distance,1,1,rotations/subarrays(ii));

error_distance = abs(estimated_distance_aux - optimal_ftm);
error_distance(~index_ap_in) = nan;

error_distance = error_distance(:);


error_distance_hf = error_distance;
error_distance_hf(~in_hf) = nan;

error_distance_lf = error_distance;
error_distance_lf(~in_lf) = nan;

subplot(3,1,2)
histogram(error_distance_hf(:), 100,'Normalization','probability')
hold on
histogram(error_distance_lf(:), 100,'Normalization','probability')
title("Distance difference")
legend("Better HF", "Better LF")
% error_all = [];

%% Disagreement
% disagreement = disagreement(:);
disagreement(disagreement == Inf) = nan;
disagreement(~index_ap_in) = nan;
disagreement = disagreement(:);
disagreement_hf = disagreement;
disagreement_hf(~in_hf) = nan;

disagreement_lf = disagreement;
disagreement_lf(~in_lf) = nan;


subplot(3,1,3),
% subplot(2,1,1)
histogram(disagreement_hf(:), 100,'Normalization','probability')
hold on
histogram(disagreement_lf(:), 100,'Normalization','probability')


legend("Better HF", "Better LF")
title("Location difference")
sgtitle("Selecting the best AP. # of arrays = " + string(subarrays(ii)))

%% Debug
% error_all = [];
% optimal_ftm_vector = optimal_ftm(:);
% dubug = [disagreement, error_angle, error_distance, index_min];
% if ( ii == 1 )
%     point_aux = 164;
%     dubug(point_aux,:)
%     optimal_ftm_vector(point_aux)
%     estimated_distance_aux(point_aux)
%     % hf
%     estimated_point_x_hf_aux = estimated_point_x_hf_rotations(:);
%     estimated_point_y_hf_aux = estimated_point_y_hf_rotations(:);
%     % lf
%     estimated_point_x_lf_aux = estimated_point_x_lf_rot(:);
%     estimated_point_y_lf_aux = estimated_point_y_lf_rot(:);
%     
%     [estimated_point_x_hf_aux(point_aux), estimated_point_y_hf_aux(point_aux)]
%     [estimated_point_x_lf_aux(point_aux), estimated_point_y_lf_aux(point_aux)]
%     disagreement_matriz = reshape(disagreement,110,5);
% 
%     hola = 1;
% end
end
end