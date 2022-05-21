clear ;
close all
clc

mkdir("../../Plots/FTM/Errors")
mkdir("../../Plots/FTM/Errors/fig")
mkdir("../../Plots/FTM/Errors/png")

% index_set_num = (1:74).';
% index_set = string(index_set_num);
% system("rm -rf mat_files/FTM/*");

routers = 4;
routers_CSI = [4,6,7,10];

load("../../mat_files/outdoor_map/data_distance_angle_true.mat")

% distances_RX = distances_RX(1:74,:);

% LOS_RX = [LOS_RX4,LOS_RX5,LOS_RX10];



load("../../mat_files/mD_track/LF/Data_3D_Direct_Path_outdoor.mat")
% angles_RX = angles_RX*(-1);

angles_dp(:,1) = angles_dp(:,1) + 18;
angles_dp(:,2) = angles_dp(:,2) + 3;
angles_dp(:,3) = angles_dp(:,3) - 3;
angles_dp(:,4) = angles_dp(:,4) - 18;

% angles_dp(1:74,3) = angles_dp(1:74,3) - 7;


error = (angles_dp - angles_RX);

% error(:,[1,2,4]) = error(:,[1,2,4]) - 1;
figure, cdfplot((error(:)))
for ii = 1:4
    hold on
    cdfplot((error(:,ii)))
end
legend("All routers", "Router 4", "Router 6", "Router 7", "Router 10", "Location", "Southeast");

figure, cdfplot(abs(error(:)))
for ii = 1:4
    hold on
    cdfplot(abs(error(:,ii)))
end
legend("All routers", "Router 4", "Router 6", "Router 7", "Router 10", "Location", "Southeast");

% error(:,5) = [];
figure, cdfplot(abs(error(:)))

% median_error = median(error,1);
% median_error([1,4]) = 2;
% median_error([2]) = 0.8;
% median_error([3]) = -0.8;
% 
% estimated_distance = estimated_distance - median_error;
% error = (estimated_distance - distances_RX);
% 
% figure, cdfplot(abs(error(:)))
% for ii = 1:routers
%     hold on
%     cdfplot(abs(error(:,ii)))
% end
% legend("All routers", "Router 4", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
% xticks(0:1:15)


% savefig("../../Plots/FTM/Errors/fig/Error");
% print("../../Plots/FTM/Errors/png/Error", "-dpng");

% 
% error_LOS = error.* LOS_RX;
% error_LOS_nans = error_LOS;
% error_LOS_nans(~LOS_RX) = nan;
% 
% error_NLOS = error.* (~LOS_RX);
% error_NLOS_nans = error_NLOS;
% error_NLOS_nans(LOS_RX) = nan;
% 
% 
% figure, cdfplot(abs(error_LOS_nans(:)))
% for ii = [1,2,5]
%     hold on
%     cdfplot(abs(error_LOS_nans(:,ii)))
% end
% legend("All routers", "Router 4", "Router 5", "Router 10", "Location", "Southeast");
% title("LoS FTM performance")
% xlabel("Error distance [m]")
% ylabel("Probability")
% 
% % savefig("../../Plots/FTM/Errors/fig/LOS_Error");
% % print("../../Plots/FTM/Errors/png/LOS_Error", "-dpng");
% 
% figure, cdfplot(abs(error_NLOS_nans(:)))
% for ii = 1:routers
%     hold on
%     cdfplot(abs(error_NLOS_nans(:,ii)))
% end
% legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
% title("NLoS FTM performance")
% xlabel("Error distance [m]")
% ylabel("Probability")

% savefig("../../Plots/FTM/Errors/fig/NLOS_Error");
% print("../../Plots/FTM/Errors/png/NLOS_Error", "-dpng");

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