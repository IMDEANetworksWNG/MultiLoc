clear
close all
clc


pwd_str = pwd;

cd ../../

mkdir("mat_files/FTM");
mkdir("mat_files/FTM/LF");


n_points = 31;
index_set_num = (1:31).';
index_set = string(index_set_num);
% system("rm -rf mat_files/FTM/*");

routers = 3;
routers_CSI = [222,223,224];


% take ground truth
% load("mat_files/Data_True/data_distance_angle_true.mat")
% load("mat_files/Data_True/Grid_LOS.mat")

% variable to save the estimated distance
estimated_distance = zeros(n_points,routers);

for id_set = 1:length(index_set_num)
    
    
    for router = 1:routers 

        load(strcat("mat_files/imdea/LF/FTM/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        estimated_distance(id_set, router) = nanmedian(FTM_distances);
    end
end

save("mat_files/FTM/LF/estimated_distance_imdea", "estimated_distance");

cd(pwd_str)