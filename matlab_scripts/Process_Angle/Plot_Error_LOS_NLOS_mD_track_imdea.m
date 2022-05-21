clear ;
close all
clc

pwd_str = pwd;

cd ../../

openfig("Map.fig")

% mkdir("Plots/AoA/Errors")
% mkdir("Plots/AoA/Errors/fig")
% mkdir("Plots/AoA/Errors/png")

index_set_num = (1:31).';
index_set = string(index_set_num);
% system("rm -rf mat_files/AoA/*");

routers = 1:3;
routers_CSI = [222;223;224];

load("mat_files/imdea_map/distance_angles.mat")
% load("mat_files/Data_True/Grid_LOS.mat")

load("mat_files/mD_track/LF/Data_3D_Direct_Path_imdea.mat")

% angles_dp = squeeze(aoa_routers(:,:,1));

angles_dp(:,3) = angles_dp(:,3) - 6;
angles_dp(:,2) = angles_dp(:,2) + 6;
angles_dp([19,20,27:31],1) = angles_dp([19,20,27:31],1)*(-1);
% 
% angles_RX = [angles_RX_4,angles_RX_5,angles_RX_10];
angles_RX = angles_RX*(-1);

error = (angles_dp - angles_RX);

figure, cdfplot(abs(error(1:18,3)))
figure, cdfplot(abs(error(2:2:12,1)))
figure, cdfplot(abs(error(1:18,1)))
figure, cdfplot(abs(error(:,2)))
% LOS_RX = [LOS_RX4,LOS_RX5,LOS_RX10];


% figure, cdfplot(abs(error(:)))
% for ii = 1:3
%     hold on
%     cdfplot(abs(error(:,ii)))
% end
% legend("All routers", "Router 4", "Router 5", "Router 10", "Location", "Southeast");
% savefig("Plots/AoA/Errors/fig/mD_track_Error");
% print("Plots/AoA/Errors/png/mD_track_Error", "-dpng");
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
% for ii = 1:3
%     hold on
%     cdfplot(abs(error_LOS_nans(:,ii)))
% end
% legend("All routers", "Router 4", "Router 5", "Router 10", "Location", "Southeast");
% title("LOS AoA performance")
% xlabel("Angle errors [deg]")
% ylabel("Probability")
% 
% savefig("Plots/AoA/Errors/fig/mD_track_LOS_Error");
% print("Plots/AoA/Errors/png/mD_track_LOS_Error", "-dpng");
% 
% figure, cdfplot(abs(error_NLOS_nans(:)))
% for ii = 1:3
%     hold on
%     cdfplot(abs(error_NLOS_nans(:,ii)))
% end
% 
% % xlim([0 25])
% legend("All routers", "Router 4", "Router 5", "Router 10", "Location", "Southeast");
% title("NLOS AoA performance")
% xlabel("Angle errors [deg]")
% ylabel("Probability")
% 
% savefig("Plots/AoA/Errors/fig/mD_track_NLOS_Error");
% print("Plots/AoA/Errors/png/mD_track_NLOS_Error", "-dpng");

cd(pwd_str)