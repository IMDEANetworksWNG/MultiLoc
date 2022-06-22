clear ;
close all
clc

index_set_num = (1:110).';
index_set = string(index_set_num);

routers = 1:5;
routers_CSI = [4,5,6,7,10];

load("mean_offset.mat")

mean_offset_final = [mean_offset(2), mean_offset(3), mean_offset(3), mean_offset(4), mean_offset(1)]
% mean_offset_final_2 = [mean_offset(2), mean_offset(3), mean_offset(3), 0, mean_offset(1)]


for id_set = 1:length(index_set_num)
    for router = 1:length(routers)  
% 
        load(strcat("../../mat_files/indoor/LF/FTM/", string(id_set), "/FTM_distances_", string(routers_CSI(router))));
        if (id_set > 74)
            mean_offset_final(4) = -(1.79 + 0.7);
%             mean_offset_final(3) = mean_offset(4);

        end
        FTM_distances = FTM_distances - mean_offset_final(router);
        save(strcat("../../mat_files/indoor/LF/FTM/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))), "FTM_distances");

    end
end