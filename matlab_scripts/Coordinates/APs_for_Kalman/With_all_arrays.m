clear 
close all
clc

%%
%pwd_str = pwd;

%cd ../../

addpath('auxiliar_functions');

load("mat_files/indoor_map/points_coordinates.mat")
load("mat_files/indoor_map/Grid_LOS.mat")
load('mat_files/indoor_map/data_distance_angle_true.mat')


load('mat_files/indoor/HF/FTM/ftm_indoor.mat')

load('mat_files/FTM/LF/estimated_distance.mat')

% load the HF data
load("mat_files/Coordinates/data_coordinates_hf_amb.mat")
% load the LF data
load("mat_files/Coordinates/data_coordinates_lf.mat")



subarrays = [8,4,2,1];
index_subarrays = {0:7,0:2:6,0:4:4,0};

estimated_distance(63,1) = 6.5;
estimated_distance(106,1) = 16.5;

estimated_point_x_hf(77,1,7) = nan;
estimated_point_y_hf(77,1,7) = nan;

router_indexes = [1 3 4 5];

calculated_distance = calculated_distance(:, router_indexes, :);
estimated_distance = estimated_distance(:, router_indexes);
estimated_point_x_hf = estimated_point_x_hf(:, router_indexes, :);
estimated_point_y_hf = estimated_point_y_hf(:, router_indexes, :);
estimated_point_x_lf = estimated_point_x_lf(:, router_indexes, :);
estimated_point_y_lf = estimated_point_y_lf(:, router_indexes, :);

[n_points, routers, rotations] = size(estimated_point_x_hf);


for ii = 1:length(subarrays)

estimated_point_x_hf_rotations = nan(n_points,routers,rotations/subarrays(ii));
estimated_point_y_hf_rotations = nan(n_points,routers,rotations/subarrays(ii));

disagreement = zeros(n_points,routers,rotations/subarrays(ii));
chosen_rotation_with_lf = nan(n_points,routers, rotations/subarrays(ii));

for id_point = 1:n_points
    for id_router = 1:routers
        for id_rotation = 1:(rotations/subarrays(ii))
            
            % take ii
            estimated_point_x_hf_aux = estimated_point_x_hf(id_point, id_router, id_rotation + index_subarrays{ii});
            estimated_point_y_hf_aux = estimated_point_y_hf(id_point, id_router, id_rotation + index_subarrays{ii});
            
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
                chosen_rotation_with_lf(id_point, id_router, id_rotation) = index_min_rotation;
                estimated_point_x_hf_rotations(id_point, id_router, id_rotation) = estimated_point_x_hf_aux;
                estimated_point_y_hf_rotations(id_point, id_router, id_rotation) = estimated_point_y_hf_aux;
            end
        end
    end
end

estimated_point_x_hf_rotations_2 = estimated_point_x_hf_rotations;
estimated_point_y_hf_rotations_2 = estimated_point_y_hf_rotations;

% save(strcat("positions_",string(subarrays(ii))), "estimated_point_x_hf_rotations_2", "estimated_point_y_hf_rotations_2");

if (ii == 4)
    th_disagreement = 3;
elseif (ii == 3)
    th_disagreement = 3.5;
else

    th_disagreement = 5;
end

estimated_point_x = zeros(n_points,routers,rotations/subarrays(ii));
estimated_point_y = zeros(n_points,routers,rotations/subarrays(ii));

for id_point = 1:n_points
    for id_router = 1:routers
        for id_rotation = 1:(rotations/subarrays(ii))
            % if the disagreement if higher that th, take LF. If not take
            % the HF data
            if (disagreement(id_point, id_router, id_rotation) > th_disagreement)
                estimated_point_x(id_point, id_router, id_rotation) = estimated_point_x_lf(id_point, id_router);
                estimated_point_y(id_point, id_router, id_rotation) = estimated_point_y_lf(id_point, id_router);
            else
                estimated_point_x(id_point, id_router, id_rotation) = estimated_point_x_hf_rotations(id_point, id_router, id_rotation);
                estimated_point_y(id_point, id_router, id_rotation) = estimated_point_y_hf_rotations(id_point, id_router, id_rotation);
            end
            % disagreement == Inf means that there is no coverage
            if (disagreement(id_point, id_router, id_rotation)  == Inf)
                estimated_point_x(id_point, id_router, id_rotation) = nan;
                estimated_point_y(id_point, id_router, id_rotation) = nan;
            end
        end
    end
end

error_lf = sqrt((points_x - estimated_point_x_lf).^2 + (points_y - estimated_point_y_lf).^2);
error_hf = sqrt((points_x - estimated_point_x_hf_rotations).^2 + (points_y - estimated_point_y_hf_rotations).^2);
error    = sqrt((points_x - estimated_point_x).^2 + (points_y - estimated_point_y).^2);

estimated_point_x_2 = estimated_point_x;
estimated_point_y_2 = estimated_point_y;

% save(strcat("positions_2_",string(subarrays(ii))), "estimated_point_x_2", "estimated_point_y_2");




% save("mat_files/Coordinates/data_coordinates_lf_hf","estimated_point_x", "estimated_point_y");
% save("mat_files/Coordinates/data_coordinates_lf_hf_rotations","estimated_point_x_hf_rotations", "estimated_point_y_hf_rotations");
% save("mat_files/Coordinates/disagreement","disagreement");

%% Plot the error of taking the AP which the HF best matches with LF

estimated_point_x_hf_lf = zeros(n_points, rotations/subarrays(ii));
estimated_point_y_hf_lf = zeros(n_points, rotations/subarrays(ii));

for id_rotation = 1:(rotations/subarrays(ii))
    
    disagreement_aux = squeeze(disagreement(:,:,id_rotation));
    estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
    % take the AP which best matches with LF
    min_disagreeement = nanmin(disagreement_aux,[],2);

    % take the index
    index_in_disagreement = disagreement_aux == min_disagreeement;

    estimated_point_x_hf_lf(:,id_rotation) = nansum(estimated_point_x_aux.*index_in_disagreement,2);
    estimated_point_y_hf_lf(:,id_rotation) = nansum(estimated_point_y_aux.*index_in_disagreement,2);

end


error_hf_hf_lf = sqrt((points_x - estimated_point_x_hf_lf).^2 + (points_y - estimated_point_y_hf_lf).^2);
% figure, cdfplot(error_hf_hf_lf(:))


%% Mean LF

estimated_point_x_lf_median = mean(estimated_point_x_lf,2);
estimated_point_y_lf_median = mean(estimated_point_y_lf,2);

error_lf_mean = sqrt((points_x - estimated_point_x_lf_median).^2 + (points_y - estimated_point_y_lf_median).^2);
figure,% cdfplot(error_lf_mean)
hold on

%% Mean LF nan

index_nan = isnan(squeeze(estimated_point_x(:,:,1)));

estimated_point_x_lf_nan = estimated_point_x_lf;
estimated_point_x_lf_nan(index_nan) = nan;

estimated_point_y_lf_nan = estimated_point_y_lf;
estimated_point_y_lf_nan(index_nan) = nan;

estimated_point_x_lf_median_nan = nanmean(estimated_point_x_lf_nan,2);
estimated_point_y_lf_median_nan = nanmean(estimated_point_y_lf_nan,2);

error_lf_mean = sqrt((points_x - estimated_point_x_lf_median_nan).^2 + (points_y - estimated_point_y_lf_median_nan).^2);
% figure, cdfplot(error_lf_mean)
% hold on

%% Distance LF

metric_distance = 1./abs(estimated_distance);
metric_distance = metric_distance./max(metric_distance,[],2);

estimated_point_x_lf_distance = sum(estimated_point_x_lf.* metric_distance,2)./sum(metric_distance,2);
estimated_point_y_lf_distance = sum(estimated_point_y_lf.* metric_distance,2)./sum(metric_distance,2);

error_lf_distance = sqrt((points_x - estimated_point_x_lf_distance).^2 + (points_y - estimated_point_y_lf_distance).^2);
hold on
%cdfplot(error_lf_distance)


%% Distance LF-2

metric = estimated_distance;
metric_2 = zeros(size(metric));

for id_point = 1:n_points
    metric_point = metric(id_point,:);
    [~, index_sort] = sort(metric_point);
    metric_point_sort = metric_point(index_sort);
    window = abs(diff(metric_point_sort([1 end])));
    metric_point_sort(index_sort);
    metric_2_aux = abs(((metric_point_sort - metric_point_sort(1)) - window)/window);
%     metric_2_aux = metric_2_aux(index_sort);

    unsorted = 1:length(metric_point);
    newInd(index_sort) = unsorted;
    
    metric_point_2 = metric_point_sort(newInd);
    hola(id_point) = isequal(metric_point, metric_point_2);
    metric_2_aux = metric_2_aux(newInd);

    metric_2(id_point,:) = metric_2_aux;
end

estimated_point_x_lf_distance_2 = sum(estimated_point_x_lf.* metric_2,2)./sum(metric_2,2);
estimated_point_y_lf_distance_2 = sum(estimated_point_y_lf.* metric_2,2)./sum(metric_2,2);

error_lf_distance_2 = sqrt((points_x - estimated_point_x_lf_distance_2).^2 + (points_y - estimated_point_y_lf_distance_2).^2);
hold on
cdfplot(error_lf_distance_2)

%save(['matlab_scripts/Coordinates/APs_for_Kalman/lf_distance_2.mat'], "estimated_point_x_lf_distance_2", "estimated_point_y_lf_distance_2");

%% Mean HF
estimated_point_x_hf_mean = nanmean(estimated_point_x_hf_rotations,2);
estimated_point_y_hf_mean = nanmean(estimated_point_y_hf_rotations,2);

error_hf_mean = sqrt((points_x - estimated_point_x_hf_mean).^2 + (points_y - estimated_point_y_hf_mean).^2);
index_nan_error_mean = isnan(error_hf_mean);
error_hf_mean(index_nan_error_mean) = Inf;
cdfplot(error_hf_mean(:))

%% Best ToF HF
% Choose the best AP and rotation based on ToF

estimated_point_x_hf_tof = nan(n_points, rotations/subarrays(ii));
estimated_point_y_hf_tof = nan(n_points, rotations/subarrays(ii));
chosen_ap_with_hf = nan(n_points, rotations/subarrays(ii));

for id_point = 1:n_points
        
    for id_rotation = 1:(rotations/subarrays(ii))
        % Select the minimum tof for AP and rotation
        candidate_subarray = nan;
        candidate_ap = nan;
        candidate_subarray_distance = Inf;

        
        for id_router = 1:routers

            distance = calculated_distance(id_point, id_router, id_rotation + index_subarrays{ii});
            
            if(~(sum(isnan(distance)) == subarrays(ii)))
                
                distance_new = ones(subarrays(1),1)*Inf;
                distance_new(id_rotation + index_subarrays{ii}) = distance;
                [dist_val, index_rotation] = nanmin(distance_new);

                if dist_val < candidate_subarray_distance

                    candidate_subarray = index_rotation;
                    candidate_ap = id_router;
                    candidate_subarray_distance = dist_val;
                    
                end
            end
        end
        % Now for this point we should have the best AP and rotation based on
        % ToF
        if(~isnan(candidate_ap) && ~isnan(candidate_subarray))
            estimated_point_x_hf_tof(id_point, id_rotation) = estimated_point_x_hf(id_point, candidate_ap, candidate_subarray);
            estimated_point_y_hf_tof(id_point, id_rotation) = estimated_point_y_hf(id_point, candidate_ap, candidate_subarray);
            chosen_ap_with_hf(id_point, id_rotation) = candidate_ap;
        end
    end
end

error_hf_tof = sqrt((points_x - estimated_point_x_hf_tof).^2 + (points_y - estimated_point_y_hf_tof).^2);
index_nan_error_tof = isnan(error_hf_tof);
error_hf_tof(index_nan_error_tof) = Inf;
cdfplot(error_hf_tof(:));

%% Plot the error of taking the AP which the HF best matches with median LF

estimated_point_x_hf_lf = zeros(n_points, rotations/subarrays(ii));
estimated_point_y_hf_lf = zeros(n_points, rotations/subarrays(ii));
chosen_ap_with_lf = nan(n_points, 1);
count = 1;

for id_rotation = 1:(rotations/subarrays(ii))
        
    
    estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
    disagreement_aux = sqrt((estimated_point_x_aux - estimated_point_x_lf_distance).^2 + (estimated_point_y_aux - estimated_point_y_lf_distance).^2);
    % take the AP which best matches with LF
    min_disagreeement = nanmin(disagreement_aux,[],2);

    % take the index
    index_in_disagreement = disagreement_aux == min_disagreeement;

    estimated_point_x_hf_lf(:,id_rotation) = nansum(estimated_point_x_aux.*index_in_disagreement,2);
    estimated_point_y_hf_lf(:,id_rotation) = nansum(estimated_point_y_aux.*index_in_disagreement,2);
    
    chosen_ap_with_lf(:, count) = sum(index_in_disagreement .* [1:routers], 2);
    count = count + 1;
end


error_hf_lf = sqrt((points_x - estimated_point_x_hf_lf).^2 + (points_y - estimated_point_y_hf_lf).^2);
%cdfplot(error_hf_lf(:))


%% Plot the error of taking the AP which the HF best matches with median LF

estimated_point_x_hf_lf_2 = zeros(n_points, rotations/subarrays(ii));
estimated_point_y_hf_lf_2 = zeros(n_points, rotations/subarrays(ii));
chosen_ap_with_lf = nan(n_points, 1);
count = 1;

for id_rotation = 1:(rotations/subarrays(ii))
        
    estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
%     estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    index_nan = isnan(squeeze(estimated_point_x(:,:,id_rotation)));
%     estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
    distances = estimated_distance;
    distances(index_nan) = nan;
    
    min_distance = nanmin(distances,[],2);
    
    index_in_disagreement = min_distance == distances;

    estimated_point_x_hf_lf_2(:,id_rotation) = nansum(estimated_point_x_aux.*index_in_disagreement,2);
    estimated_point_y_hf_lf_2(:,id_rotation) = nansum(estimated_point_y_aux.*index_in_disagreement,2);
    
    index_0 = estimated_point_x_hf_lf_2(:,id_rotation) == 0 | estimated_point_y_hf_lf_2(:,id_rotation) == 0;
    if (sum(index_0) > 0)
        
        estimated_point_x_hf_lf_2(index_0, id_rotation) = estimated_point_x_lf_distance_2(index_0);
        estimated_point_y_hf_lf_2(index_0, id_rotation) = estimated_point_y_lf_distance_2(index_0);
        
        estimated_point_x(index_0,1,id_rotation) = estimated_point_x_lf_distance_2(index_0);
        estimated_point_y(index_0,1,id_rotation) = estimated_point_y_lf_distance_2(index_0);
    end
%     estimated_point_x_hf_lf_2(:,id_rotation)
    
    chosen_ap_with_lf(:, count) = sum(index_in_disagreement .* [1:routers], 2);
    count = count + 1;
end


error_hf_lf_2 = sqrt((points_x - estimated_point_x_hf_lf_2).^2 + (points_y - estimated_point_y_hf_lf_2).^2);

%% check if other estimates are closer, so that we can average it

estimated_point_x_hf_lf_new = zeros(size(estimated_point_x_hf_lf));
estimated_point_y_hf_lf_new = zeros(size(estimated_point_x_hf_lf));

th = 1;

for id_point = 1:n_points

    for id_rotation = 1:(rotations/subarrays(ii))

        estimated_point_x_aux = squeeze(estimated_point_x(id_point,:,id_rotation));
        estimated_point_y_aux = squeeze(estimated_point_y(id_point,:,id_rotation));
        
        distance_estimates = sqrt( (estimated_point_x_hf_lf_2(id_point,id_rotation) - estimated_point_x_aux).^2 + (estimated_point_y_hf_lf_2(id_point,id_rotation) - estimated_point_y_aux).^2);
        index_in = distance_estimates <= th;
        
        estimated_point_x_hf_lf_new(id_point, id_rotation) = nanmean(estimated_point_x_aux(index_in));
        estimated_point_y_hf_lf_new(id_point, id_rotation) = nanmean(estimated_point_y_aux(index_in));

    end
end
% 
error_hf_hf_lf_new = sqrt((points_x - estimated_point_x_hf_lf_new).^2 + (points_y - estimated_point_y_hf_lf_new).^2);
save(['matlab_scripts/Coordinates/APs_for_Kalman/for_kalman_' num2str(subarrays(ii)) '.mat'], "estimated_point_x_hf_lf_new", "estimated_point_y_hf_lf_new");


cdfplot(error_hf_lf_2(:))
cdfplot(error_hf_hf_lf_new(:))

% cdfplot(error_hf_lf_2(:))


% save('with_8_arrays.mat', 'error_hf_lf', 'error_lf_mean', 'error_hf_mean')

% Add the optimal
load('mat_files/indoor/HF/optimal_rotation_index.mat');

% Only LOS
load('mat_files/indoor_map/Grid_LOS.mat')

estimated_point_optimal_x_hf = nan(n_points, routers);
estimated_point_optimal_y_hf = nan(n_points, routers);

for ap_id=1:routers
    for point_id=1:n_points
        
        if LOS_RX(point_id, ap_id) == 1
            estimated_point_optimal_x_hf(point_id, ap_id) = estimated_point_x_hf(point_id, ap_id, indoor_optimal_rotation_index(point_id, ap_id));
            estimated_point_optimal_y_hf(point_id, ap_id) = estimated_point_y_hf(point_id, ap_id, indoor_optimal_rotation_index(point_id, ap_id));
        end
    end
end

error_hf_optimal = sqrt((points_x - estimated_point_optimal_x_hf).^2 + (points_y - estimated_point_optimal_y_hf).^2);

% cdfplot(error_hf_hf_lf(:))


% cdfplot(error_hf_optimal(:))
title(strcat("# of subarrays ",string(subarrays(ii)))); 
legend("LF-dist2","HF mean", "HF ToF", "LF-HF", "LF-HF-New")

if subarrays(ii) == 1
    xlim([0 26])
    xticks(1:26)

elseif subarrays(ii) == 2
    xlim([0 18])
    xticks(1:18)

elseif subarrays(ii) == 4
    xlim([0 16])
    xticks(1:16)
elseif subarrays(ii) == 8
    xlim([0 16])
    xticks(1:16)
end

%save(['matlab_scripts/Coordinates/Only_2_APs/all_indoor_' num2str(subarrays(ii)) '.mat'], "error_lf_mean", "error_lf_distance", "error_lf_distance_2", "error_hf_mean", "error_hf_lf", "error_hf_hf_lf_new", "error_hf_tof");

%% Are we choosing the same in LF and HF?
% chosen_ap_with_lf
% chosen_ap_with_hf
load('mat_files/indoor_map/data_distance_angle_true.mat', 'AP_labels_HF')


% figure
% plot(1:110, chosen_ap_with_lf, 'o', 'LineStyle', 'none');
% hold on
% plot(1:110, chosen_ap_with_hf, 'x', 'LineStyle', 'none');
% 
% yticks([1:5])
% yticklabels(string(AP_labels_HF))
nans = isnan(error_hf_hf_lf_new);
sum(sum(nans))
end

% %% Further testing
% 
% % Load FTM
% load('mat_files/indoor/HF/FTM/ftm_indoor.mat')
% 
% 
% % Chose between LF and HF but take that one which has minimum ToF
% % FTM LF
% load('mat_files/FTM/LF/estimated_distance.mat')
% 
% chosen_rotation = nan(n_points, routers);
% 
% % Here we want the rotation that has lower ToF
% for id_ap=1:routers
%     for id_point=1:n_points
%         
%         % Get LF ToF
%         lf_tof = estimated_distance(id_point, id_ap);
%         hf_tof = nan(subarrays(ii), 1);
% 
%         for id_rotation=1:subarrays(ii)
% 
%             hf_tof = calculated_distance(id_point, id_ap,id_rotation + index_subarrays{ii})
%         end
% %         
% %         for id_rotation=1:subarrays(ii)
% %         
% %             hf_tof(id_rotation) = calculated_distance(id_point, id_ap, id_rotation);
% %         end
%         
%         [min_hf, index] = nanmin(hf_tof);
%         
%         %if min_hf <= lf_tof
%         chosen_rotation(id_point, id_ap) = index;
%         %end
%     end
% end
% 
% chosen_ap = nan(110, 1);
% 
% % Now we need to chose the best AP for HF
% for id_point=1:110
%     
%     hf_tof = nan(5, 1);
%     
%     for id_ap=1:5
%         
%         hf_tof(id_ap) = calculated_distance(id_point, id_ap, chosen_rotation(id_point, id_ap));
%     end
%     
%     [min_hf, index] = nanmin(hf_tof);
%     
%     chosen_ap(id_point) = index;
% end
% 
% 
% % estimated_point_x_lf_median = mean(estimated_point_x_lf,2);
% % estimated_point_y_lf_median = mean(estimated_point_y_lf,2);
% 
% % LF
% tof_lf_hf_x = mean(estimated_point_x_lf,2);
% tof_lf_hf_y = mean(estimated_point_y_lf,2);
% 
% % HF
% tof_lf_hf_x = nan(110, 1);
% tof_lf_hf_y = nan(110, 1);
% 
% % Uncomment the following line if you prefer to chose an AP based on LF
% %chosen_ap = chosen_ap_with_lf;
% 
% for id_point=1:110
%         
%     tof_lf_hf_x(id_point) = estimated_point_x_hf(id_point, chosen_ap(id_point), chosen_rotation(id_point, chosen_ap(id_point)));
%     tof_lf_hf_y(id_point) = estimated_point_y_hf(id_point, chosen_ap(id_point), chosen_rotation(id_point, chosen_ap(id_point)));
% end
% 
% error_hf_lf_tof = sqrt((points_x - tof_lf_hf_x).^2 + (points_y - tof_lf_hf_y).^2);
% 
% cdfplot(error_hf_lf_tof(:))
% 
% % HF-LF
% hf_lf_difference = sqrt((estimated_point_x_lf_median - tof_lf_hf_x).^2 + (estimated_point_y_lf_median - tof_lf_hf_y).^2);
% 
% diff_threshold = 5;
% 
% lf_hf_estimated_x = nan(110, 1);
% lf_hf_estimated_y = nan(110, 1);
% 
% for id_point=1:110
%     
%     % Take LF if difference too big
%     if hf_lf_difference(id_point) > diff_threshold
%         
%         lf_hf_estimated_x(id_point) = estimated_point_x_lf_median(id_point);
%         lf_hf_estimated_y(id_point) = estimated_point_y_lf_median(id_point);
% 
%     else %Take 
%         
%         lf_hf_estimated_x(id_point) = tof_lf_hf_x(id_point);
%         lf_hf_estimated_y(id_point) = tof_lf_hf_y(id_point);      
%     end
% end
% 
% error_hf_lf_tof = sqrt((points_x - lf_hf_estimated_x).^2 + (points_y - lf_hf_estimated_y).^2);
% 
% h = cdfplot(error_hf_lf_tof(:));
% set(h, 'LineStyle', '--');
% 
% h = cdfplot(error_hf_lf_tof(:));
% set(h, 'LineStyle', '-', 'LineWidth',2);
% 
% legend("Mean LF", "Mean HF", "LF-HF", "LF-HF2", "HF optimal rotation", "HF based on ToF", "New LF/HF", "All")
% title("Estimated position by all APs, 8 arrays")
% % end