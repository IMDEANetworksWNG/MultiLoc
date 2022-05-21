clear ;
close all
clc

pwd_str = pwd;

cd ../../

mkdir("Plots/")
mkdir("Plots/AoA/")
mkdir("Plots/AoA/fig")
mkdir("Plots/AoA/png")

index_set_num = (1:20).';
index_set = string(index_set_num);
% system("rm -rf mat_files/AoA/*");

routers = 2;
routers_CSI = [4,5,6,7,10];

load("mat_files/in_out_map/data_distance_angle_true.mat")
% load("mat_files/indoor_map/Grid_LOS.mat")

load("mat_files/in_out/HF/CSI/csi_in_out2.mat")
load("mat_files/in_out/HF/FTM/ftm_in_out.mat")

% load("mat_files/index_optimal_rotation")


% angles_dp = aoa_optimal_rotation;


% angles_RX = angles_RX*(-1);


[n_points, routers, rotationss] = size(calculated_aoa);
% angles_RX(15,:) = [];

index_optimal_rotation = zeros(n_points, routers);
for ii = 1:n_points
    for jj = 1:routers
       [~, index_optimal_rotation(ii,jj)] = nanmin(calculated_distance(ii,jj,:)); 
    end
end

angles_dp = zeros(n_points, routers);

for ii = 1:n_points
    for jj = 1:routers
        angles_dp(ii,jj) = calculated_aoa(ii,jj,index_optimal_rotation(ii,jj));
    end
end
angles_RX = angles_RX(1:n_points,:)*(-1);
error = angles_dp - angles_RX;


% figure%, cdfplot(abs(error(:)))
% for ii = 1:routers
%     hold on
%     cdfplot((error(:,ii)))
% end

points = 1:10;
figure, cdfplot(abs(error(points,1)+1))

points = 11:20;

figure, cdfplot(abs(error(points,2)+2.2))
% legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
% savefig("Plots/AoA/fig/mD_track_Error");
% print("Plots/AoA/png/mD_track_Error", "-dpng");

% error_LOS = error.* LOS_RX;
% error_LOS_nans = error_LOS;
% error_LOS_nans(~LOS_RX) = nan;
% 
% error_NLOS = error.* (~LOS_RX);
% error_NLOS_nans = error_NLOS;
% error_NLOS_nans(LOS_RX) = nan;
% 
% 
% % error_LOS_nans = error_LOS_nans - aoa_offset;
% figure, cdfplot(abs(error_LOS_nans(:)))
% for ii = [1,2,5]
%     hold on
%     cdfplot(abs(error_LOS_nans(:,ii)))
% end
% legend("All routers", "Router 4", "Router 5", "Router 10", "Location", "Southeast");title("LOS AoA performance")
% xlabel("Angle errors [deg]")
% ylabel("Probability")

% savefig("Plots/AoA/fig/mD_track_LOS_Error");
% print("Plots/AoA/png/mD_track_LOS_Error", "-dpng");

% figure, cdfplot(abs(error_NLOS_nans(:)))
% for ii = 1:routers
%     hold on
%     cdfplot(abs(error_NLOS_nans(:,ii)))
% end
% 
% % xlim([0 25])
% legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");title("NLOS AoA performance")
% xlabel("Angle errors [deg]")
% ylabel("Probability")

% savefig("Plots/AoA/fig/mD_track_NLOS_Error");
% print("Plots/AoA/png/mD_track_NLOS_Error", "-dpng");

% aoa_offset = nanmedian(error_LOS_nans,1);
% aoa_offset(3) = median(error_NLOS_nans([32:36 107:108],3));
% aoa_offset(4) = median(error_NLOS_nans([70:74 109:110],4));
% save("mat_files/indoor/HF/CSI/aoa_offset", "aoa_offset");


cd(pwd_str)