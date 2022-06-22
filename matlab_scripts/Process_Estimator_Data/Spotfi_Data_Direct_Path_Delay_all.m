clear
clc
close all

pwd_str = pwd;

cd ../../
mkdir("mat_files/")
mkdir("mat_files/Jade_all")

scenarios = ["indoor", "outdoor", "in_out","imdea"];

routers_all = {1:5,1:4,1:2,1:3};
index_set_num_all = {(1:110).',(111:126).',(1:20).',(1:31).'};
for scenario = 1:length(scenarios)
strcat("mat_files/Jade_all/Data_2D_5_120_processed_Delay_",scenarios(scenario),".mat")
    load(strcat("mat_files/Jade_all/Data_2D_5_120_processed_Delay_",scenarios(scenario),".mat"))

    routers = routers_all{scenario};
    index_set_num = index_set_num_all{scenario};
    angles_dp = zeros(length(index_set_num),length(routers));


    for id_point = 1:length(index_set_num)
        for id_router = 1:length(routers)

            angles_dp(id_point,id_router) = aoa_routers{id_point,id_router}(1,1);

        end 
    end


    save(strcat("mat_files/Jade_all/Data_2D_5_120_Direct_Path_Delay_",scenarios(scenario),".mat"), "angles_dp")

end
cd(pwd_str)