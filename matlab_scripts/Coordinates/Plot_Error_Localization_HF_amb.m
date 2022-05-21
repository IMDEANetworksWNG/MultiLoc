clear 
close all
clc
% load the coordinates

% pwd_str = pwd;

% cd ../../

mkdir("Plots/Coordinates/")
mkdir("Plots/Coordinates/fig")
mkdir("Plots/Coordinates/png")

load("mat_files/indoor_map/Grid_LOS.mat")


% load the ground truth
load("mat_files/indoor_map/points_coordinates.mat")

% load the estimated coordinates
load("mat_files/Coordinates/data_coordinates_hf_amb")
[n_points, routers, rotations] = size(estimated_point_x_hf);



%% Take the optimal
% load the index to the optimal 
load("mat_files/indoor/HF/optimal_rotation_index.mat")

estimated_point_x_optimal = zeros(n_points,routers);
estimated_point_y_optimal = zeros(n_points,routers);

for id_point = 1:n_points
    for id_router = 1:routers
        index_optimal = indoor_optimal_rotation_index(id_point,id_router);
        estimated_point_x_optimal(id_point, id_router) = estimated_point_x_hf(id_point, id_router,index_optimal);
        estimated_point_y_optimal(id_point, id_router) = estimated_point_y_hf(id_point, id_router,index_optimal);
        
    end
end

% take the error
error = sqrt((points_x - estimated_point_x_optimal).^2 + (points_y - estimated_point_y_optimal).^2);

figure, cdfplot(abs(error(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
% savefig("Plots/Coordinates/fig/Error");
% print("Plots/Coordinates/png/Error", "-dpng");


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
cdfplot(abs(error_NLOS_nans([32:36 107:108],3)))
cdfplot(abs(error_NLOS_nans([70:74 109:110],4)))

legend("All routers", "Router 4", "Router 5", "Router 10", "Router 6", "Router 7", "Location", "Southeast");
title("LoS Coordinates performance")
xlabel("Error distance [m]")
ylabel("Probability")

% savefig("Plots/Coordinates/fig/LOS_Error");
% print("Plots/Coordinates/png/LOS_Error", "-dpng");

figure, cdfplot(abs(error_NLOS_nans(:)))
for ii = 1:routers
    hold on
    cdfplot(abs(error_NLOS_nans(:,ii)))
end
legend("All routers", "Router 4", "Router 5", "Router 6", "Router 7", "Router 10", "Location", "Southeast");
title("NLoS Coordinates performance")
xlabel("Error distance [m]")
ylabel("Probability")

% savefig("Plots/Coordinates/fig/NLOS_Error");
% print("Plots/Coordinates/png/NLOS_Error", "-dpng");

% cd(pwd_str)