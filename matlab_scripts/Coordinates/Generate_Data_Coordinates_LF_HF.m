clear 
close all
clc

pwd_str = pwd;

cd ../../

load("mat_files/indoor_map/points_coordinates.mat")
load("mat_files/indoor_map/Grid_LOS.mat")


% load the HF data
load("mat_files/Coordinates/data_coordinates_hf.mat")
% load the LF data
load("mat_files/Coordinates/data_coordinates_lf.mat")

[n_points, routers, rotations] = size(estimated_point_x_hf);

estimated_point_x_hf_rotations = nan(n_points,routers,rotations/2);
estimated_point_y_hf_rotations = nan(n_points,routers,rotations/2);

disagreement = zeros(n_points,routers,rotations/2);
chosen_rotation_with_lf = nan(n_points,routers, rotations/2);

for id_point = 1:n_points
    for id_router = 1:routers
        for id_rotation = 1:(rotations/2)
            
            % first subarray
            estimated_point_x_hf_aux(1) = estimated_point_x_hf(id_point, id_router, id_rotation);
            estimated_point_y_hf_aux(1) = estimated_point_y_hf(id_point, id_router, id_rotation);

            % second subarray
            estimated_point_x_hf_aux(2) = estimated_point_x_hf(id_point, id_router, id_rotation + 4);
            estimated_point_y_hf_aux(2) = estimated_point_y_hf(id_point, id_router, id_rotation + 4);
            
            % take the lf one
            estimated_point_x_lf_aux = estimated_point_x_lf(id_point, id_router);
            estimated_point_y_lf_aux = estimated_point_y_lf(id_point, id_router);
            
            % take the one that better agree with LF
            distance_arrays = sqrt((estimated_point_x_hf_aux - estimated_point_x_lf_aux).^2 + (estimated_point_y_hf_aux - estimated_point_y_lf_aux).^2);
            
            % check if we have coverage, if not take the lf information
            if (sum(isnan(distance_arrays)) == 2)
                disagreement(id_point, id_router, id_rotation) = Inf;
            else
                % if we have coverage, take the subarray that better agrees
                % with LF
%                 [~, index_min_rotation] = min(distance_arrays);
                [~, index_min_rotation] = min(squeeze(calculated_distance(id_point, id_router,:)));

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

% figure, cdfplot(disagreement(1:25,5,1))
% hold on
% cdfplot(disagreement(37:61,2,1))
% figure
% cdfplot(disagreement(1:25,5,2))
% hold on
% cdfplot(disagreement(37:61,2,2))
% figure
% cdfplot(disagreement(1:25,5,3))
% hold on
% cdfplot(disagreement(37:61,2,3))
% figure
% cdfplot(disagreement(1:25,5,4))
% hold on
% cdfplot(disagreement(37:61,2,4))

th_disagreement = 2;

estimated_point_x = zeros(n_points,routers,rotations/2);
estimated_point_y = zeros(n_points,routers,rotations/2);

for id_point = 1:n_points
    for id_router = 1:routers
        for id_rotation = 1:(rotations/2)
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

save("mat_files/Coordinates/data_coordinates_lf_hf","estimated_point_x", "estimated_point_y");
save("mat_files/Coordinates/data_coordinates_lf_hf_rotations","estimated_point_x_hf_rotations", "estimated_point_y_hf_rotations");
save("mat_files/Coordinates/disagreement","disagreement");

% figure, 
% cdfplot(error_lf(1:25,5))
% hold on
% cdfplot(error_hf(1:25,5,:))
% cdfplot(error(1:25,5,:))
% legend("LF", "HF", "LF-HF")
% 
% figure, 
% cdfplot(error_lf(37:61,2))
% hold on
% cdfplot(error_hf(37:61,2,:))
% cdfplot(error(37:61,2,:))
% legend("LF", "HF", "LF-HF")

%% Plot the error of taking the AP which the HF best matches with LF

estimated_point_x_hf_lf = zeros(n_points, rotations/2);
estimated_point_y_hf_lf = zeros(n_points, rotations/2);

for id_rotation = 1:(rotations/2)
    
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
figure, cdfplot(error_lf_mean)
hold on

%% Mean HF
estimated_point_x_hf_mean = nanmean(estimated_point_x_hf_rotations,2);
estimated_point_y_hf_mean = nanmean(estimated_point_y_hf_rotations,2);

error_hf_mean = sqrt((points_x - estimated_point_x_hf_mean).^2 + (points_y - estimated_point_y_hf_mean).^2);

cdfplot(error_hf_mean(:))

%% Plot the error of taking the AP which the HF best matches with median LF

estimated_point_x_hf_lf = zeros(n_points, rotations/2);
estimated_point_y_hf_lf = zeros(n_points, rotations/2);


for id_rotation = 1:(rotations/2)
        
    
    estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
    disagreement_aux = sqrt((estimated_point_x_aux - estimated_point_x_lf_median).^2 + (estimated_point_y_aux - estimated_point_y_lf_median).^2);
    % take the AP which best matches with LF
    min_disagreeement = nanmin(disagreement_aux,[],2);

    % take the index
    index_in_disagreement = disagreement_aux == min_disagreeement;

    estimated_point_x_hf_lf(:,id_rotation) = nansum(estimated_point_x_aux.*index_in_disagreement,2);
    estimated_point_y_hf_lf(:,id_rotation) = nansum(estimated_point_y_aux.*index_in_disagreement,2);

end


error_hf_lf = sqrt((points_x - estimated_point_x_hf_lf).^2 + (points_y - estimated_point_y_hf_lf).^2);
cdfplot(error_hf_lf(:))

legend("Mean LF", "Mean HF", "LF-HF")
title("Estimated position by all APs")


%% Plot the error of taking the AP which the HF best matches with median LF

estimated_point_x_hf_lf = zeros(n_points, rotations/2);
estimated_point_y_hf_lf = zeros(n_points, rotations/2);


for id_rotation = 1:(rotations/2)
        
    
    estimated_point_x_aux = squeeze(estimated_point_x_hf(:,:,id_rotation));
    estimated_point_y_aux = squeeze(estimated_point_y_hf(:,:,id_rotation));
    disagreement_aux = sqrt((estimated_point_x_aux - estimated_point_x_lf_median).^2 + (estimated_point_y_aux - estimated_point_y_lf_median).^2);
    % take the AP which best matches with LF
    min_disagreeement = nanmin(disagreement_aux,[],2);

    estimated_point_x_aux = squeeze(estimated_point_x(:,:,id_rotation));
    estimated_point_y_aux = squeeze(estimated_point_y(:,:,id_rotation));
    
    % take the index
    index_in_disagreement = disagreement_aux == min_disagreeement;

    estimated_point_x_hf_lf(:,id_rotation) = nansum(estimated_point_x_aux.*index_in_disagreement,2);
    estimated_point_y_hf_lf(:,id_rotation) = nansum(estimated_point_y_aux.*index_in_disagreement,2);

end


error_hf_lf = sqrt((points_x - estimated_point_x_hf_lf).^2 + (points_y - estimated_point_y_hf_lf).^2);
figure, 
cdfplot(error_hf_lf(:))

% %% combine LOS + NLOS
% 
% estimated_point_x_all = squeeze(mean(estimated_point_x,2, "omitnan"));
% estimated_point_y_all = squeeze(mean(estimated_point_y,2,"omitnan"));
% 
% error_all = sqrt((points_x - estimated_point_x_all).^2 + (points_y - estimated_point_y_all).^2);
% 
% figure, cdfplot(error_all(:))
% 
% 
cd(pwd_str)

