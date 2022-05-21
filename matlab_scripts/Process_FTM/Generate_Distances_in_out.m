clear
close all
clc


pwd_str = pwd;

cd ../../

mkdir("mat_files/FTM");
mkdir("mat_files/FTM/LF");


n_points = 20;
index_set_num = (1:n_points).';
index_set = string(index_set_num);
% system("rm -rf mat_files/FTM/*");

routers = 2;
routers_CSI = [10,4];

% variable to save the estimated distance
estimated_distance = zeros(n_points, routers);

for id_set = 1:length(index_set_num)
    
    
    for router = 1:routers 

        load(strcat("mat_files/in_out/LF/FTM/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        estimated_distance(id_set, router) = median(FTM_distances);
    end
end

save("mat_files/FTM/LF/estimated_distance_in_out", "estimated_distance");

cd(pwd_str)