close all
clear
clc

pwd_str = pwd;

cd ../../

for n_e = 3

load(strcat("mat_files/Jade/Data_1D_AoA_",string(n_e),".mat"))
% load("mat_files/data_distance_angle_true.mat")

% mkdir("mat_files/Active_Localization")

routers = [4;5;6;7;10];
index_set_num = (1:110).';

angles_dp = zeros(length(index_set_num), length(routers));


for ii = 1:length(index_set_num)
    for jj = 1:length(routers)
        
        aoa_aux = AoA{ii,jj};
        
        if (isempty(aoa_aux))
            aoa_aux = nan;
            power_aux = nan;
            
        else
            aoa_aux = AoA{ii,jj}{1,1};
            power_aux = power{ii,jj}{1,1};
            [~, index_max] = max(power_aux);
            angles_dp(ii,jj) = aoa_aux(index_max);
        end
    end 
end



save(strcat("mat_files/Jade/Data_MUSIC_AoA_Direct_Path_",string(n_e)), "angles_dp");    
end
cd(pwd_str)