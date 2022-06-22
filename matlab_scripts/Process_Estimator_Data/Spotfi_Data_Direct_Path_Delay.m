clear
clc
close all

pwd_str = pwd;

cd ../../

load("mat_files/Jade/Data_2D_5_120_processed_Delay.mat")
% load("mat_files/data_distance_angle_true.mat")

% mkdir("mat_files/Active_Localization")

routers = [4;5;6;7;10];
index_set_num = (1:110).';

angles_dp = zeros(74,length(routers));


for id_point = 1:length(index_set_num)
    for id_router = 1:length(routers)

        angles_dp(id_point,id_router) = aoa_routers{id_point,id_router}(1,1);

    end 
end



save("mat_files/Jade/Data_2D_5_120_Direct_Path_Delay", "angles_dp");    

cd(pwd_str)