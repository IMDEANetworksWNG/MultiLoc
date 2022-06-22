clear
clc
close all

pwd_str = pwd;

cd ../../

load("mat_files/Jade/Data_2D_5_120_processed_Delay_outdoor2.mat")
% load("mat_files/data_distance_angle_true.mat")

% mkdir("mat_files/Active_Localization")

index_set_num = (111:126).';
index_set = string(index_set_num);

routers = [4;6;7;10];

angles_dp = zeros(length(index_set_num),length(routers));


for id_point = 1:length(index_set_num)
    for id_router = 1:length(routers)

        angles_dp(id_point,id_router) = aoa_routers{id_point,id_router}(1,1);

    end 
end



save("mat_files/Jade/Data_2D_5_120_Direct_Path_Delay_outdoor2", "angles_dp");    

cd(pwd_str)