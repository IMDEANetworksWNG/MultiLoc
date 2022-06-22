close all
clear
clc

pwd_str = pwd;

cd ../../

load("mat_files/mD_track/LF/Data_3D_processed_in_out.mat")
% load("mat_files/data_distance_angle_true.mat")

% mkdir("mat_files/Active_Localization")

routers = [10;4];
index_set_num = (1:20).';

angles_dp = zeros(20,length(routers));


for id_point = 1:length(index_set_num)
    for id_router = 1:length(routers)
        angles_dp(id_point,id_router) = aoa_routers{id_point,id_router}(1,1);
    end 
end



save("mat_files/mD_track/LF/Data_3D_Direct_Path_in_out", "angles_dp");    

cd(pwd_str)