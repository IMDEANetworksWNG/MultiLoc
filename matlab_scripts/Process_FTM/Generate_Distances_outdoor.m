clear
close all
clc


pwd_str = pwd;

cd ../../

% mkdir("mat_files/FTM");


index_set_num = (111:126).';
index_set = string(index_set_num);
% system("rm -rf mat_files/FTM/*");

routers = 4;
routers_CSI = [4,6,7,10];


% take ground truth
% load("mat_files/Data_True/Grid_LOS.mat")

% variable to save the estimated distance
estimated_distance = zeros(length(index_set), routers);

for id_set = 1:length(index_set_num)
    
    
    for router = 1:routers 

        load(strcat("mat_files/outdoor/LF/FTM/", string(index_set_num(id_set)), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        estimated_distance(id_set, router) = median(FTM_distances);
    end
end

save("mat_files/FTM/LF/estimated_distance_outdoor", "estimated_distance");

cd(pwd_str)