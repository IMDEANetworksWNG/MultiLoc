clear ;
close all
clc

index_set_num = (111:126).';
index_set = string(index_set_num);

routers = 1:4;
routers_CSI = [4,6,7,10];

load("mean_offset.mat")

mean_offset_final = [mean_offset(2), mean_offset(4) + 0.3, -(1.79 + 0.7), (mean_offset(1) + 0.7)]
% mean_offset_final_2 = [mean_offset(2), mean_offset(3), mean_offset(3), 0, mean_offset(1)]


for id_set = 1:length(index_set_num)
    for router = 1:length(routers)  
        strcat("../../mat_files/outdoor/LF/FTM/", string(index_set_num(id_set)), "/FTM_distances_", string(routers_CSI(router)))
        load(strcat("../../mat_files/outdoor/LF/FTM/", string(index_set_num(id_set)), "/FTM_distances_", string(routers_CSI(router))));
        FTM_distances = FTM_distances - mean_offset_final(router);
        save(strcat("../../mat_files/outdoor/LF/FTM/", string(index_set_num(id_set)), "/FTM_distances_calibrated_", string(routers_CSI(router))), "FTM_distances");

    end
end