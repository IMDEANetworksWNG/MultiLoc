clear ;
close all
clc

mkdir("../../Plots/FTM/")
mkdir("../../Plots/FTM/fig")
mkdir("../../Plots/FTM/png")

index_set_num = (1:110).';
index_set = string(index_set_num);
% system("rm -rf mat_files/FTM/*");

routers = 5;
routers_CSI = [4,5,6,7,10];

load("../../mat_files/indoor_map/data_distance_angle_true.mat")
load("../../mat_files/indoor_map/Grid_LOS.mat")

% LOS_RX = [LOS_RX4,LOS_RX5,LOS_RX10];

% distances_RX = [distances_RX_4,distances_RX_5,distances_RX_10];
% estimated_distance = zeros(size(distances_RX));


load("../../mat_files/indoor/HF/FTM/ftm_optimal_rotation")

error = (ftm_optimal_rotation - distances_RX);

figure, cdfplot(abs(error(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
% savefig("../../Plots/FTM/fig/Error");
% print("../../Plots/FTM/png/Error", "-dpng");


error_LOS = error.* LOS_RX;
error_LOS_nans = error_LOS;
error_LOS_nans(~LOS_RX) = nan;

error_NLOS = error.* (~LOS_RX);
error_NLOS_nans = error_NLOS;
error_NLOS_nans(LOS_RX) = nan;


figure, cdfplot(abs(error_LOS_nans(:)))
for ii = [1,2,5]
    hold on
    cdfplot(abs(error_LOS_nans(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
title("LoS FTM performance")
xlabel("Error distance [m]")
ylabel("Probability")

% savefig("../../Plots/FTM/fig/LOS_Error");
% print("../../Plots/FTM/png/LOS_Error", "-dpng");

figure, cdfplot(abs(error_NLOS_nans(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error_NLOS_nans(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
title("NLoS FTM performance")
xlabel("Error distance [m]")
ylabel("Probability")

% savefig("../../Plots/FTM/fig/NLOS_Error");
% print("../../Plots/FTM/png/NLOS_Error", "-dpng");

% 
% load("../../mat_files/Data_True/Grid_NLOS.mat")
% 
% NLOS_1_RX = [NLOS_1_RX4,NLOS_1_RX5,NLOS_1_RX6,NLOS_1_RX10];
% NLOS_2_RX = [NLOS_2_RX4,NLOS_2_RX5,NLOS_2_RX6,NLOS_2_RX10];
% 
% 
% 
% error_NLOS_1 = nan(size(error));
% error_NLOS_1(NLOS_1_RX) = error(NLOS_1_RX);
% 
% error_NLOS_2 = nan(size(error));
% error_NLOS_2(NLOS_2_RX) = error(NLOS_2_RX);
% 
% 
% figure, cdfplot(error_NLOS_1(:))
% for ii = 1:4
%     hold on
%     cdfplot(error_NLOS_1(:,ii))
% end
% legend("All routers", "Router 4", "Router 5", "Router 6", "Router 10", "Location", "Southeast");
% title("1 wall NLOS FTM performance")
% xlabel("Error distance [m]")
% ylabel("Probability")
% 
% savefig("../../Plots/FTM/Errors/fig/NLOS_1_wall_Error");
% print("../../Plots/FTM/Errors/png/NLOS_1_wall_Error", "-dpng");
% 
% figure, cdfplot(error_NLOS_2(:))
% for ii = 1:4
%     hold on
%     cdfplot(error_NLOS_2(:,ii))
% end
% legend("All routers", "Router 4", "Router 5", "Router 6", "Router 10", "Location", "Southeast");
% title("More than 1 wall NLOS FTM performance")
% xlabel("Error distance [m]")
% ylabel("Probability")
% 
% savefig("../../Plots/FTM/Errors/fig/NLOS_2_wall_Error");
% print("../../Plots/FTM/Errors/png/NLOS_2_wall_Error", "-dpng");
% 
% figure, plot(error_NLOS_1)
% legend(string([4;5;6;10]))
% 
% figure, plot(error_NLOS_2)
% legend(string([4;5;6;10]))