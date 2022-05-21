clear 
close all
clc
% load the coordinates

pwd_str = pwd;

cd ../../

mkdir("Plots/Coordinates/")
mkdir("Plots/Coordinates/fig")
mkdir("Plots/Coordinates/png")

load("mat_files/indoor_map/Grid_LOS.mat")

routers = 5;

% load the ground truth
load("mat_files/indoor_map/points_coordinates.mat")

% load the estimated coordinates
load("mat_files/Coordinates/LF/data_coordinates")

% take the error
error = sqrt((points_x - estimated_point_x).^2 + (points_y - estimated_point_y).^2)

figure, cdfplot(abs(error(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
savefig("Plots/Coordinates/fig/Error");
print("Plots/Coordinates/png/Error", "-dpng");


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
title("LoS Coordinates performance")
xlabel("Error distance [m]")
ylabel("Probability")

savefig("Plots/Coordinates/fig/LOS_Error");
print("Plots/Coordinates/png/LOS_Error", "-dpng");

figure, cdfplot(abs(error_NLOS_nans(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error_NLOS_nans(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
title("NLoS Coordinates performance")
xlabel("Error distance [m]")
ylabel("Probability")

savefig("Plots/Coordinates/fig/NLOS_Error");
print("Plots/Coordinates/png/NLOS_Error", "-dpng");

cd(pwd_str)