close all
clear
clc

pwd_str = pwd;

cd ../../

load("mat_files/Jade/Data_2D_5_120_outdoor.mat")
% load("mat_files/data_distance_angle_true.mat")

% mkdir("mat_files/Active_Localization")

% index_set = 300;
index_set_num = (111:126).';
index_set = string(index_set_num);

routers = [4;6;7;10];

direct_angle_mean = zeros(length(index_set_num), length(routers));
direct_angle_std = zeros(length(index_set_num), length(routers));
direct_angle = {};


aoa_routers = cell(size(AoA));
toa_routers = cell(size(AoA));
power_routers = cell(size(AoA));

for ii = 1:length(index_set_num)
    for jj = 1:length(routers)
        aoa_aux = AoA{ii,jj};
%         aod_aux = AoD{ii,jj};
        power_aux = power{ii,jj};

        toa_aux = ToA{ii,jj};
        aoa_dp = [];
        for kk = 1

            aoa_ind = aoa_aux{kk,1}(1:end);
            power_ind = power_aux{kk,1}(1:end);
            toa_ind = toa_aux{kk,1}(1:end);
            toa_ind = toa_ind - min(toa_ind);

            [~, index_sort] = sort(toa_ind, "ascend");
            power_ind = power_ind(index_sort);
            aoa_ind = aoa_ind(index_sort);
            toa_ind = toa_ind(index_sort);
            if (length(toa_ind) > 1)

                [~, index_sort] = sort(power_ind(1:2), "descend");

                
                if (isequal(index_sort, [2 1].') && diff(toa_ind(1:2)) > 5)
                    index_sort = [1 2];
                    in = 1
                end
                
                aoa_ind(1:2) = aoa_ind(index_sort);
                power_ind(1:2) = power_ind(index_sort);
                toa_ind(1:2) = toa_ind(index_sort);
                
                aoa_dp = aoa_ind(index_sort(1));
                toa_dp = toa_ind(index_sort(1));
                power_dp = power_ind(index_sort(1));
                aod_dp = 0;
            else
                aoa_ind = aoa_ind(1);
                toa_ind = toa_ind(1);
                power_ind = power_ind(1);
                aod_dp = 0;
            end

        end
        aoa_routers{ii,jj} = aoa_ind;
        toa_routers{ii,jj} = toa_ind;
        power_routers{ii,jj} = power_ind;
        

    end 
end



save("mat_files/Jade/Data_2D_5_120_processed_Delay_outdoor2", "aoa_routers", "toa_routers","power_routers");    

cd(pwd_str)