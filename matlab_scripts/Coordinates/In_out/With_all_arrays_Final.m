clear 
close all
clc

%%
pwd_str = pwd;

cd ../../../

addpath('auxiliar_functions');

load("mat_files/indoor_map/points_coordinates.mat")
load("mat_files/indoor_map/Grid_LOS.mat")

load('mat_files/indoor/HF/FTM/ftm_indoor.mat')

load('mat_files/FTM/LF/estimated_distance.mat')

% load the HF data
load("mat_files/Coordinates/data_coordinates_hf_amb.mat")
% load the LF data
load("mat_files/Coordinates/data_coordinates_lf.mat")

[n_points, routers, rotations] = size(estimated_point_x_hf);


subarrays = [8,4,2,1];
index_subarrays = {0:7,0:2:6,0:4:4,0};
% 
% estimated_distance(63,1) = 6.5;
% estimated_distance(106,1) = 16.5;

estimated_point_x_hf(77,1,7) = nan;
estimated_point_y_hf(77,1,7) = nan;

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

figure,
hold on

%% Distance LF

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

error_lf_distance = sqrt((points_x - estimated_point_x_lf_distance_2).^2 + (points_y - estimated_point_y_lf_distance_2).^2);
hold on
cdfplot(error_lf_distance)

% save(['matlab_scripts/Coordinates/With_multiple_antenna_arrays/lf_distance_2.mat'], "estimated_point_x_lf_distance_2", "estimated_point_y_lf_distance_2");


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
cdfplot(error_hf_tof(:));



%% MultiLoc

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


cdfplot(error_hf_lf_2(:))
% cdfplot(error_hf_hf_lf_new(:))

error_hf_hf_lf_new = error_hf_lf_2;
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
legend("LF-dist","HF mean", "HF ToF", "LF-HF")

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

save(['matlab_scripts//Coordinates/With_multiple_antenna_arrays/all_indoor_' num2str(subarrays(ii)) '.mat'], "error_lf_distance", "error_hf_hf_lf_new", "error_hf_tof");

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
cd(pwd_str)
