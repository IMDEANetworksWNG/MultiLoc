clear 
close all
clc

pwd_str = pwd;

cd ../../

load("mat_files/indoor_map/points_coordinates.mat")
load("mat_files/indoor_map/Grid_LOS.mat")

load('mat_files/indoor/HF/FTM/ftm_indoor.mat')


% load the HF data
load("mat_files/Coordinates/data_coordinates_hf.mat")
% load the index of the best rotation for HF data
load("mat_files/indoor/HF/optimal_rotation_index.mat")
% load the LF data
load("mat_files/Coordinates/data_coordinates_lf.mat")
% load the LF/HF data
load("mat_files/Coordinates/data_coordinates_lf_hf")


[n_points, routers, rotations] = size(estimated_point_x_hf);

%% Plot best LF
estimated_point_x_lf_best = sum(estimated_point_x_lf .* LOS_RX,2);
estimated_point_y_lf_best = sum(estimated_point_y_lf .* LOS_RX,2);

% points in the small room on the left
index_in = [32:36 107:108];
router_in = 3;
estimated_point_x_lf_best(index_in) = estimated_point_x_lf(index_in,router_in);
estimated_point_y_lf_best(index_in) = estimated_point_y_lf(index_in,router_in);

% points in the small room on the right
index_in = [70:74 109:110];
router_in = 4;
estimated_point_x_lf_best(index_in) = estimated_point_x_lf(index_in,router_in);
estimated_point_y_lf_best(index_in) = estimated_point_y_lf(index_in,router_in);

error_lf_best = sqrt((points_x - estimated_point_x_lf_best).^2 + (points_y - estimated_point_y_lf_best).^2);
figure, cdfplot(error_lf_best)
title("Best AP. LF")

%% Plot mean LF
estimated_point_x_lf_mean = mean(estimated_point_x_lf,2);
estimated_point_y_lf_mean = mean(estimated_point_y_lf,2);



error_lf_mean = sqrt((points_x - estimated_point_x_lf_mean).^2 + (points_y - estimated_point_y_lf_mean).^2);
figure, cdfplot(error_lf_mean)
title("Mean APs. LF")
%% Plot median LF
estimated_point_x_lf_median = median(estimated_point_x_lf,2);
estimated_point_y_lf_median = median(estimated_point_y_lf,2);



error_lf_median = sqrt((points_x - estimated_point_x_lf_median).^2 + (points_y - estimated_point_y_lf_median).^2);
figure, cdfplot(error_lf_median)
title("Median APs. LF")
%% Plot all LF
% estimated_point_x_lf_mean = mean(estimated_point_x_lf,2);
% estimated_point_y_lf_mean = mean(estimated_point_y_lf,2);



error_lf = sqrt((points_x - estimated_point_x_lf).^2 + (points_y - estimated_point_y_lf).^2);
figure, cdfplot(error_lf(:))
title("All APs. LF")

%% Plot best rotation and best AP hF
estimated_point_x_hf_aux = zeros(size(estimated_point_x_lf));
estimated_point_y_hf_aux = zeros(size(estimated_point_x_lf));
for id_point = 1:n_points
    for id_router = 1:routers
        estimated_point_x_hf_aux(id_point,id_router) = estimated_point_x_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
        estimated_point_y_hf_aux(id_point,id_router) = estimated_point_y_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
    end
end

estimated_point_x_hf_best = nansum(estimated_point_x_hf_aux .* LOS_RX,2);
estimated_point_y_hf_best = nansum(estimated_point_y_hf_aux .* LOS_RX,2);

% points in the small room on the left
index_in = [32:36 107:108];
router_in = 3;
estimated_point_x_hf_best(index_in) = estimated_point_x_hf_aux(index_in,router_in);
estimated_point_y_hf_best(index_in) = estimated_point_y_hf_aux(index_in,router_in);

% points in the small room on the right
index_in = [70:74 109:110];
router_in = 4;
estimated_point_x_hf_best(index_in) = estimated_point_x_hf_aux(index_in,router_in);
estimated_point_y_hf_best(index_in) = estimated_point_y_hf_aux(index_in,router_in);

error_hf_best = sqrt((points_x - estimated_point_x_hf_best).^2 + (points_y - estimated_point_y_hf_best).^2);
figure, cdfplot(error_hf_best)
title("Best rotation and best AP. HF")

%% Plot best rotation and mean AP hF
estimated_point_x_hf_aux = zeros(size(estimated_point_x_lf));
estimated_point_y_hf_aux = zeros(size(estimated_point_x_lf));
for id_point = 1:n_points
    for id_router = 1:routers
        estimated_point_x_hf_aux(id_point,id_router) = estimated_point_x_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
        estimated_point_y_hf_aux(id_point,id_router) = estimated_point_y_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
    end
end

estimated_point_x_hf_best_mean = mean(estimated_point_x_hf_aux,2, "omitnan");
estimated_point_y_hf_best_mean = mean(estimated_point_y_hf_aux,2, "omitnan");

error_hf_best_mean = sqrt((points_x - estimated_point_x_hf_best_mean).^2 + (points_y - estimated_point_y_hf_best_mean).^2);
figure, cdfplot(error_hf_best_mean)
title("Best rotation and mean APs. HF")

%% Plot all rotations all aps

error_hf_hf = sqrt((points_x - estimated_point_x_hf).^2 + (points_y - estimated_point_y_hf).^2);
figure, cdfplot(error_hf_hf(:))
title("All rotations and all APs. HF")

% %% Plot only 4 rotations considering 2 subarrays all aps
% 
% load("mat_files/Coordinates/data_coordinates_lf_hf_rotations");
% 
% 
% error_hf_lf_4_rotations = sqrt((points_x - estimated_point_x_hf_rotations).^2 + (points_y - estimated_point_y_hf_rotations).^2);
% figure, cdfplot(error_hf_lf_4_rotations(:))
% 
% %% Plot best rotation and AP wighted mean std of the power
% 
% load("mat_files/LOS_NLOS/std_all")
% 
% 
% estimated_point_x_hf_aux = zeros(size(estimated_point_x_lf));
% estimated_point_y_hf_aux = zeros(size(estimated_point_x_lf));
% std_all_aux = zeros(size(estimated_point_x_lf));
% 
% for id_point = 1:n_points
%     for id_router = 1:routers
%         estimated_point_x_hf_aux(id_point,id_router) = estimated_point_x_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
%         estimated_point_y_hf_aux(id_point,id_router) = estimated_point_y_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
%         std_all_aux(id_point,id_router) = std_all(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
%     end
% end
% 
% std_all_aux = 1./std_all_aux;
% 
% estimated_point_x_hf_best_weigthed_std = nansum(estimated_point_x_hf_aux .* std_all_aux,2)./nansum(std_all_aux,2);
% estimated_point_y_hf_best_weigthed_std = nansum(estimated_point_y_hf_aux .* std_all_aux,2)./nansum(std_all_aux,2);
% 
% 
% error_hf_best_weigthed_std = sqrt((points_x - estimated_point_x_hf_best_weigthed_std).^2 + (points_y - estimated_point_y_hf_best_weigthed_std).^2);
% figure, cdfplot(error_hf_best_weigthed_std)
% 
% %% Plot best rotation and AP which better matches with LF
% 
% estimated_point_x_hf_aux = zeros(size(estimated_point_x_lf));
% estimated_point_y_hf_aux = zeros(size(estimated_point_x_lf));
% disagreement = zeros(size(estimated_point_x_lf));
% 
% for id_point = 1:n_points
%     for id_router = 1:routers
%         estimated_point_x_hf_aux(id_point,id_router) = estimated_point_x_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
%         estimated_point_y_hf_aux(id_point,id_router) = estimated_point_y_hf(id_point,id_router, indoor_optimal_rotation_index(id_point, id_router));
%         disagreement(id_point,id_router) = sqrt((estimated_point_x_hf_aux(id_point,id_router) - estimated_point_x_lf(id_point,id_router)).^2 + (estimated_point_y_hf_aux(id_point,id_router) - estimated_point_y_lf(id_point,id_router)).^2);
%     end
% end
% 
% % take the AP which best matches with LF
% min_disagreeement = nanmin(disagreement,[],2);
% 
% % take the index
% index_in_disagreement = disagreement == min_disagreeement;
% 
% estimated_point_x_hf_best_lf = nansum(estimated_point_x_hf_aux.*index_in_disagreement,2);
% estimated_point_y_hf_best_lf = nansum(estimated_point_y_hf_aux.*index_in_disagreement,2);
% 
% error_hf_best_best_lf = sqrt((points_x - estimated_point_x_hf_best_lf).^2 + (points_y - estimated_point_y_hf_best_lf).^2);
% figure, cdfplot(error_hf_best_best_lf)
% 
% %% Plot best rotation and AP which better matches with LF for the four rotations
% 
% load("mat_files/Coordinates/data_coordinates_lf_hf","estimated_point_x", "estimated_point_y");
% load("mat_files/Coordinates/disagreement","disagreement");
% 
% estimated_point_x_hf_lf = zeros(n_points, rotations/2);
% estimated_point_y_hf_lf = zeros(n_points, rotations/2);
% 
% 
% for id_rotation = 1:(rotations/2)
%     
%     disagreement_aux = squeeze(disagreement(:,:,id_rotation));
%     estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
%     estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
%     % take the AP which best matches with LF
%     min_disagreeement = nanmin(disagreement_aux,[],2);
% 
%     % take the index
%     index_in_disagreement = disagreement_aux == min_disagreeement;
% 
%     estimated_point_x_hf_lf(:,id_rotation) = nansum(estimated_point_x_aux.*index_in_disagreement,2);
%     estimated_point_y_hf_lf(:,id_rotation) = nansum(estimated_point_y_aux.*index_in_disagreement,2);
% 
% end
% 
% 
% error_hf_hf_lf = sqrt((points_x - estimated_point_x_hf_lf).^2 + (points_y - estimated_point_y_hf_lf).^2);
% figure, cdfplot(error_hf_hf_lf(:))

cd(pwd_str)