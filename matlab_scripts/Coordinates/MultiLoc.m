clear 
% close all
clc

%%
pwd_str = pwd;

cd ../../

mkdir("plots")
mkdir("plots/MultiLoc")

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
mkdir("mat_files/MultiLoc")
        
th_disagreement_all = [4.5,4.5,4.5,4.5];


for scenario = 1:5
% for scenario = 2
mkdir(strcat("mat_files/MultiLoc/", scenarios_save(scenario)))
% load the coordinates of the client
load(strcat("mat_files/",scenarios(scenario),"_map/points_coordinates.mat"))

% load the HF FTM data
load(strcat("mat_files/",scenarios(scenario),"/HF/FTM/ftm_",scenarios(scenario),".mat"))
load(strcat("mat_files/",scenarios(scenario),"/HF/CSI/csi_",scenarios(scenario),".mat"))
calculated_aoa = -1*calculated_aoa;
if (scenario == 1 || scenario == 3 )
    load(strcat("mat_files/",scenarios(scenario),"/HF/calculated_aoa_amb.mat"))
    calculated_aoa = angles_dp;
end

if (scenario == 5)
    load(strcat("mat_files/",scenarios(scenario),"/HF/CSI/angles_dp_amb.mat"))
    angles_dp(:,1,:) = angles_dp(:,1,:) + 3;
    angles_dp(:,3,:) = angles_dp(:,3,:) - 14;
    calculated_aoa = angles_dp;
end

% load the LF FTM data
load(strcat("mat_files/FTM/LF/estimated_distance_",scenarios(scenario),".mat"))
load(strcat("mat_files/mD_track/LF/Data_3D_Direct_Path_",scenarios(scenario),".mat"))


% load the HF data
load(strcat("mat_files/Coordinates/data_coordinates_hf_",scenarios(scenario),".mat"))
% load the LF data
load(strcat("mat_files/Coordinates/data_coordinates_lf_",scenarios(scenario),".mat"))


[n_points, routers, rotations] = size(estimated_point_x_hf);

fig = figure;

if (scenario == 2)
    points_x = points_x.';
    points_y = points_y.';
    estimated_distance(:, [1 4]) = estimated_distance(:, [4 1]);

elseif (scenario == 3)
    router_indexes = [2 5];
    calculated_distance = calculated_distance(:, router_indexes, :);
    estimated_distance = estimated_distance(:, router_indexes);
    angles_dp = angles_dp(:, router_indexes);
    calculated_aoa = calculated_aoa(:, router_indexes, :);
    estimated_point_x_hf = estimated_point_x_hf(:, router_indexes, :);
    estimated_point_y_hf = estimated_point_y_hf(:, router_indexes, :);
    estimated_point_x_lf = estimated_point_x_lf(:, router_indexes, :);
    estimated_point_y_lf = estimated_point_y_lf(:, router_indexes, :);

    [n_points, routers, rotations] = size(estimated_point_x_hf);
    
end
alpha = 0.85;
for ii = 1:length(subarrays)

estimated_point_x_hf_rotations = nan(n_points,routers,rotations/subarrays(ii));
estimated_point_y_hf_rotations = nan(n_points,routers,rotations/subarrays(ii));

disagreement = zeros(n_points,routers,rotations/subarrays(ii));
calculated_distance_rotation = nan(n_points,routers,rotations/subarrays(ii));
calculated_aoa_rotation = nan(n_points,routers,rotations/subarrays(ii));

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
        
            % take the lf one
            calculated_distance_aux = calculated_distance(id_point, id_router, id_rotation + index_subarrays{ii});
            calculated_aoa_aux = calculated_aoa(id_point, id_router, id_rotation + index_subarrays{ii});

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
                
                calculated_distance_rotation(id_point, id_router, id_rotation) = calculated_distance_aux(index_min_rotation);
                calculated_aoa_rotation(id_point, id_router, id_rotation) = calculated_aoa_aux(index_min_rotation);

            end
        end
    end
end

estimated_point_x_hf_rotations_2 = estimated_point_x_hf_rotations;
estimated_point_y_hf_rotations_2 = estimated_point_y_hf_rotations;

% save(strcat("positions_",string(subarrays(ii))), "estimated_point_x_hf_rotations_2", "estimated_point_y_hf_rotations_2");
th_disagreement = th_disagreement_all(ii);
% if (ii == 4)
%     th_disagreement = 3;
% elseif (ii == 3)
%     th_disagreement = 3.5;
% else
% 
%     th_disagreement = 4.5;
% end

estimated_point_x = zeros(n_points,routers,rotations/subarrays(ii));
estimated_point_y = zeros(n_points,routers,rotations/subarrays(ii));

for id_point = 1:n_points
    for id_router = 1:routers
        for id_rotation = 1:(rotations/subarrays(ii))
            % if the disagreement if higher that th, take LF. If not take
            % the HF data
            rule_mmWave_larger = (estimated_distance(id_point, id_router) - calculated_distance_rotation(id_point, id_router, id_rotation)) < -1.5;
            rule_aoa = abs(angles_dp(id_point, id_router) - calculated_aoa_rotation(id_point, id_router, id_rotation)) > 0;
            
            rule_1 = rule_mmWave_larger && rule_aoa;
            rule_1 = rule_mmWave_larger;
            
            rule_disagreement = disagreement(id_point, id_router, id_rotation) > th_disagreement;
            takes_lf = rule_disagreement || rule_1;
%             takes_lf = rule_disagreement;
            if (takes_lf)
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

estimated_distance_rotations = repmat(estimated_distance,1,1,rotations/subarrays(ii));
% save(strcat("positions_2_",string(subarrays(ii))), "estimated_point_x_2", "estimated_point_y_2");




% save("mat_files/Coordinates/data_coordinates_lf_hf","estimated_point_x", "estimated_point_y");
% save("mat_files/Coordinates/data_coordinates_lf_hf_rotations","estimated_point_x_hf_rotations", "estimated_point_y_hf_rotations");
% save("mat_files/Coordinates/disagreement","disagreement");
subplot(2,2,ii)
hold on

%% MultiLoc LF

metric = estimated_distance;
metric_2 = zeros(size(metric));
newInd = [];

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

error_multiloc_lf = sqrt((points_x - estimated_point_x_lf_distance_2).^2 + (points_y - estimated_point_y_lf_distance_2).^2);
hold on
h = cdfplot2(error_multiloc_lf(:));%, "colors", colors(1,:))
set(h, "color",  colors(1,:))

%% MultiLoc ToF
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

error_multiloc_hf = sqrt((points_x - estimated_point_x_hf_tof).^2 + (points_y - estimated_point_y_hf_tof).^2);
h = cdfplot2(error_multiloc_hf(:));%, "colors", colors(1,:))
set(h, "color",  colors(2,:))


%% MultiLoc

estimated_point_x_hf_lf_2 = zeros(n_points, rotations/subarrays(ii));
estimated_point_y_hf_lf_2 = zeros(n_points, rotations/subarrays(ii));
chosen_ap_with_lf = nan(n_points, 1);
count = 1;

index_ap_in = zeros( n_points, routers, rotations/subarrays(ii));

for id_rotation = 1:(rotations/subarrays(ii))
        
    estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
%     estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    index_nan = isnan(squeeze(estimated_point_x(:,:,id_rotation)));
%     estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
    distance_metric = alpha*estimated_distance;
    angles_dp_aux = abs(angles_dp);
    angles_dp_aux(angles_dp_aux < 45) = 0;
    angle_metric = (1-alpha)*angles_dp_aux;
    distances = distance_metric + angle_metric;
    index_below = angles_dp_aux < 45;
    distances(index_below) = estimated_distance(index_below);
    distances(index_nan) = nan;
    
    min_distance = nanmin(distances,[],2);
    index_in_disagreement = min_distance == distances;
    index_ap_in(:,:,id_rotation) = index_in_disagreement;
    % check if matches with HF distance
%     for id_point = 1:n_points
%         hf_distance = calculated_distance(id_point,:,id_rotation);
%         diff_aux = abs(min_distance(id_point,1) - hf_distance);
%         if (sum(diff_aux < 5) > 0)
%             
%             actual_ap = find(index_in_disagreement(id_point,:) == 1);
%             [~,next_ap] = min(diff_aux);
% %             next_ap = find(diff_aux > 5,1 ,'first');
%             index_ap_in(id_point,actual_ap,id_rotation) = false;
%             index_ap_in(id_point,next_ap,id_rotation) = true;
%         end
%         
%     end
    

    
    estimated_point_x_hf_lf_2(:,id_rotation) = nansum(estimated_point_x_aux.*index_in_disagreement,2);
    estimated_point_y_hf_lf_2(:,id_rotation) = nansum(estimated_point_y_aux.*index_in_disagreement,2);
    
    index_0 = estimated_point_x_hf_lf_2(:,id_rotation) == 0 & estimated_point_y_hf_lf_2(:,id_rotation) == 0;
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


error_multiloc = sqrt((points_x - estimated_point_x_hf_lf_2).^2 + (points_y - estimated_point_y_hf_lf_2).^2);


% cdfplot2(error_multiloc(:))
h = cdfplot2(error_multiloc(:));%, "colors", colors(1,:))
set(h, "color",  colors(3,:))



title(strcat("# of subarrays ",string(subarrays(ii)))); 
legend("MultiLoc LF","MultiLoc HF", "MultiLoc")


sgtitle(scenarios_save(scenario))
%% Saving the data
save(strcat("mat_files/MultiLoc/",scenarios_save(scenario),"/all_",num2str(subarrays(ii)),".mat"), "error_multiloc_lf", "error_multiloc_hf", "error_multiloc");
savefig(fig,strcat("plots/MultiLoc/",scenarios_save(scenario)))
end
end


cd(pwd_str)