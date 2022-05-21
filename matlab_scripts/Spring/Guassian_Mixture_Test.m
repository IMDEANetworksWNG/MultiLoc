clear 
close all
clc

pwd_str = pwd;

cd ../../

load("mat_files/Spring/paths.mat")

[n_points, routers] = size(paths);
routers_CSI = [4,5,6,7,10];

estimated_distance = zeros(size(paths));
    
% load the data
for id_set = 20
    
    
    for router = 5 

        load(strcat("mat_files/FTM/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        % take the estimated distance as the mean
        index_neg = FTM_distances < 0;
        FTM_distances(index_neg) = median(FTM_distances);
        estimated_distance_point = log(FTM_distances);
        
        % take the # of paths
        paths_point = paths(id_set,router);
        % fit the data into a guassian mixture model
        try
            GMModel = fitgmdist(estimated_distance_point,paths_point);
        catch
            GMModel = fitgmdist(estimated_distance_point,1);
        end
        % take the mean of the first guassian
        mus = GMModel.mu;
        estimated_distance(id_set, router) = exp(min(mus));
        
        
    end
end

save("mat_files/Spring/estimated_distance_spring", "estimated_distance")

cd(pwd_str)