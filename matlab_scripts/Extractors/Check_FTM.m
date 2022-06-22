clear ;
close all
clc

index_set_num = (1:74).';
index_set = string(index_set_num);

routers = 1:5;
routers_CSI = [4,5,6,7,10];

load("mean_offset.mat")

mean_offset_final = [mean_offset(2), mean_offset(3), mean_offset(3), mean_offset(4), mean_offset(1)];

arethesame = false(74,5);

for id_set = 1:length(index_set_num)
    for router = 1:length(routers)  
        load(strcat("../../mat_files/FTM/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        FTM_distances_1 = FTM_distances;
        
        load(strcat("../../mat_files/FTM_2/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        FTM_distances_2 = FTM_distances;
        arethesame(id_set,router) = isequal(FTM_distances_1,FTM_distances_2);

    end
end

ok = sum(arethesame(:))/numel(arethesame)