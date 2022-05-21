clear
close all
clc


pwd_str = pwd;

cd ../../

mkdir("mat_files/FTM");
mkdir("mat_files/FTM/LF");


n_points = 110;
index_set_num = (1:n_points).';
index_set = string(index_set_num);
% system("rm -rf mat_files/FTM/*");

routers = 5;
routers_CSI = [4,5,6,7,10];

% variable to save the estimated distance
estimated_distance = zeros(n_points, routers);

for id_set = 1:length(index_set_num)
    
    
    for router = 1:routers 

        load(strcat("mat_files/indoor/LF/FTM/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        estimated_distance(id_set, router) = mean(FTM_distances);
    end
end

save("mat_files/FTM/LF/estimated_distance_indoor", "estimated_distance");

cd(pwd_str)