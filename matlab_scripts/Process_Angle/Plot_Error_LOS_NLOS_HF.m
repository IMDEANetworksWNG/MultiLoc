clear ;
close all
clc

% pwd_str = pwd;

% cd ../../

mkdir("Plots/")
mkdir("Plots/AoA/")
mkdir("Plots/AoA/fig")
mkdir("Plots/AoA/png")

index_set_num = (1:110).';
index_set = string(index_set_num);
% system("rm -rf mat_files/AoA/*");

routers = 5;
routers_CSI = [4,5,6,7,10];

load("mat_files/indoor_map/data_distance_angle_true.mat")
load("mat_files/indoor_map/Grid_LOS.mat")

load("mat_files/indoor/HF/CSI/aoa_optimal_rotation.mat")

angles_dp = aoa_optimal_rotation;


% angles_RX = angles_RX*(-1);

error = angles_dp - angles_RX;


figure, cdfplot(abs(error(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
% savefig("Plots/AoA/fig/mD_track_Error");
% print("Plots/AoA/png/mD_track_Error", "-dpng");

error_LOS = error.* LOS_RX;
error_LOS_nans = error_LOS;
error_LOS_nans(~LOS_RX) = nan;

error_NLOS = error.* (~LOS_RX);
error_NLOS_nans = error_NLOS;
error_NLOS_nans(LOS_RX) = nan;


% error_LOS_nans = error_LOS_nans - aoa_offset;
figure, cdfplot(abs(error_LOS_nans(:)))
for ii = [1,2,5]
    hold on
    cdfplot(abs(error_LOS_nans(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 10", "Location", "Southeast");title("LOS AoA performance")
xlabel("Angle errors [deg]")
ylabel("Probability")

% savefig("Plots/AoA/fig/mD_track_LOS_Error");
% print("Plots/AoA/png/mD_track_LOS_Error", "-dpng");

figure, cdfplot(abs(error_NLOS_nans(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error_NLOS_nans(:,ii)))
end

% xlim([0 25])
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");title("NLOS AoA performance")
xlabel("Angle errors [deg]")
ylabel("Probability")

% savefig("Plots/AoA/fig/mD_track_NLOS_Error");
% print("Plots/AoA/png/mD_track_NLOS_Error", "-dpng");

aoa_offset = nanmedian(error_LOS_nans,1);
aoa_offset(3) = median(error_NLOS_nans([32:36 107:108],3));
aoa_offset(4) = median(error_NLOS_nans([70:74 109:110],4));
save("mat_files/indoor/HF/CSI/aoa_offset", "aoa_offset");


% cd(pwd_str)