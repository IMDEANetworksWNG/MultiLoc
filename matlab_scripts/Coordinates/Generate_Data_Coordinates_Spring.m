clear 
% close all
clc

pwd_str = pwd;

cd ../../

% load the MUSIC data
load("mat_files/Jade/Data_MUSIC_AoA_Direct_Path_3.mat")
angles_dp(:,2) = angles_dp(:,2) - 3.5;
angles_dp(26:38,5) = angles_dp(26:38,5) - 5;
angles_dp(:,3) = angles_dp(:,3) - 7;

% load the FTM ranging
load("mat_files/Spring/estimated_distance_spring.mat")

% load the true data
% load("mat_files/Data_True/points_coordinates.mat")
load("mat_files/indoor_map/RX_coordinates.mat")
EEso 
[n_points, routers] = size(angles_dp);

% mkdir("Plots/Paper/Localization")
% mkdir("Plots/Paper/Localization/Spring/")
% mkdir("Plots/Paper/Localization/Spring/png/")
% mkdir("Plots/Paper/Localization/Spring/fig/")
% mkdir("Plots/Paper/Localization/Spring/pdf/")

%% Estimate the positions
for id_router = 1:routers
    for id_point = 1:n_points

        estimated_point_x_p(id_point,id_router) = sind(angles_dp(id_point,id_router))*estimated_distance(id_point,id_router);
        estimated_point_y_p(id_point,id_router) = cosd(angles_dp(id_point,id_router))*estimated_distance(id_point,id_router);


%         figure, plot(estimated_point_x_p(id_point,id_router), estimated_point_y_p(id_point,id_router), "*")
        

        if (id_router == 1 || id_router == 3 || id_router == 4)
            estimated_point_x_spring(id_point,id_router) = RX_x(id_router) - (estimated_point_x_p(id_point,id_router));
            estimated_point_y_spring(id_point,id_router) = RX_y(id_router) - (estimated_point_y_p(id_point,id_router)*(-1));
        else
            estimated_point_x_spring(id_point,id_router) = RX_x(id_router) - (estimated_point_x_p(id_point,id_router)*(-1));
            estimated_point_y_spring(id_point,id_router) = RX_y(id_router) - estimated_point_y_p(id_point,id_router);
        end
    end
end

save("mat_files/Coordinates/data_coordinates_spring","estimated_point_x_spring", "estimated_point_y_spring");

cd(pwd_str)