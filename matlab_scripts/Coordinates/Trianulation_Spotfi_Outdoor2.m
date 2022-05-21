clear 
clc
% close all



pwd_str = pwd;

cd ../../

addpath("functions")

load("mat_files/outdoor_map/data_distance_angle_true.mat")
load("mat_files/outdoor_map/points_coordinates.mat")
% load("mat_files/walls.mat")
load("mat_files/outdoor_map/RX_coordinates.mat")
openfig("mat_files/outdoor_map/outdoor_map")

load("mat_files/Jade/Data_2D_5_120_Direct_Path_Delay_outdoor2.mat")

% angles_dp = squeeze(aoa_routers(:,:,1));
% angles_dp = angles_RX;

angles_dp(:,[1,4]) = angles_dp(:,[4,1]);

% angles_dp(:,1) = angles_dp(:,1) + 18;
angles_dp(:,2) = angles_dp(:,2) + 3;
angles_dp(:,3) = angles_dp(:,3) - 3;
% angles_dp(:,4) = angles_dp(:,4) - 11;

% Calibrate
% angles_dp(:,2) = angles_dp(:,2) + 3.5;
% angles_dp(:,3) = angles_dp(:,3) - 2;

% mkdir("Plots/Paper/Localization")
% mkdir("Plots/Paper/Localization/triangulation/")
% mkdir("Plots/Paper/Localization/triangulation/png/")
% mkdir("Plots/Paper/Localization/triangulation/fig/")
% mkdir("Plots/Paper/Localization/triangulation/pdf/")

% merge them
% RX_x = [RX_x_4; RX_x_5; RX_x_10];
% RX_y = [RX_y_4; RX_y_5; RX_y_10];

% aps_position = [RX_x, RX_y];

aps_position = [RX_x.', RX_y.'];
% plot them
% figure, plot(aps(:,1),aps(:,2),"*")
% 
% number of points
[n_points, ~] = size(angles_RX);

% move to [x,y] coordinates
points = [points_x.', points_y.'];

% routers_CSI = [4, 5, 5,7 10];
routers = 4;

% put the max and min number to create the grid
min_angle = [-90 -90 -90 -90 -90];
max_angle = [90 90 90 90 90];

% inverse the symbol
% angles_dp = (-1)*angles_dp;
% angles_dp(:,4) = (-1)*angles_dp(:,4);
% angles_dp(:,2) = (-1)*angles_dp(:,2);

% Plot_Map_No_Points()

limit_x = 17.9;
limit_y = 24.4;


rotation_routers = [135, -135, -45, 180];

for id_point = 1:n_points
    
    % take the position to test
    position = points(id_point,:);
%     plot(position(1),position(2),"*");
    
    position_total(id_point,:) = position;
    % and plot it
%     hold on
%     plot(position(1), position(2), "*")

    % estimate each angle
    angles_point = angles_dp(id_point,:);
    
    
    % heatmap for ap 1
    % vector for x and Y coordinates
    grid_vector_x = 0:0.1:limit_x;
    grid_vector_y = 0:0.1:limit_y;
    
    % matrices X and Y
    [X,Y] = meshgrid(grid_vector_x,grid_vector_y);

    % heatmap for all aps
    heatmap_total = nan([size(X),3]);
    
    % loop to iterate over aps
    for id_router = 1:routers


        % take the angles to test
        angles_test = min_angle(id_router):0.1:max_angle(id_router);

        % heatmap for the specific ap
        heatmap = zeros(size(X));

        % standar deviation for the Gaussian function
        std_gauss = 3;
        % calculate the power distrubution in terms of the angle of the
        % direct path
        power_dist = exp(-((angles_test-angles_point(id_router)).^2)/(2*(std_gauss^2)));
        
        % compute the th to iterate over all the angles to test as many of
        % them has really low power
        th = 0.05;
        index_remove = power_dist <= th;
        
        % remove the position with the lowest power values
        power_dist(index_remove) = [];
        angles_test(index_remove) = [];
        
        % estimate the grid where each position has an angle
        diff_X = (X - RX_x(id_router));
        diff_Y = (Y - RX_y(id_router));
        
        % take the first decimal
        angles_ap = round(rad2deg(angle(diff_X + 1i*diff_Y)),1);    
        angles_ap = angles_ap - rotation_routers(id_router);
        angles_ap = round(angles_ap,1);
%         figure, pcolor(angles_ap)
%         colorbar

        
        angles_test = round(angles_test,1);
%         if (id_router == 1)
%         figure, pcolor(angles_ap)
%         colorbar
%         angles_ap = angles_ap + 45;
%         figure, pcolor(angles_ap)
%         colorbar
%             angles_test = angles_test + 45;
%         end
        
%         if (id_router == 2)
%             figure, pcolor(angles_ap)
%             colorbar
%             angles_ap = angles_ap +135;
%             figure, pcolor(angles_ap)
%             colorbar
%             angles_test = angles_test + 45;
%         end
        
        % loop to assing to each angle to searcht its power
        isinside = zeros(length(angles_test),1);
        for id_angle = 1:length(angles_test)
            
            index_change = angles_ap == angles_test(id_angle);
            isinside(id_angle,1) = sum(index_change(:));
            heatmap(index_change) = power_dist(id_angle);
     
        end
        
        % if the heatmap is 0 means that there is no power, but it is
        % better to set as the th
        index_0 = heatmap == 0;
        heatmap(index_0) = th;
        
        figure, pcolor(heatmap)
        colorbar

        % save to a big matriz
        heatmap_total(:,:,id_router) = heatmap;
    end


    % take the product of the heatmaps
    heatmap = prod(heatmap_total,3);
    
%     figure, pcolor(heatmap)
%     colorbar
%     
    % take the maximum
    [max_num,max_idx] = max(heatmap(:));  
    % find it in the matrix
    [row, column] = find(heatmap == max_num);
    % take the position
    solution = [X(row(1), column(1)),Y(row(1), column(1))];
    solution_total(id_point,:) = solution;
    
    % estimate the error
    error(id_point,1) = sqrt(sum(((solution) - position).^2,2));
%     close all
end
[~, index_max] = max(error);
% error(index_max) = 10;
fig_error = figure;
cdfplot(error);
save("mat_files/CDF_data/outdoor/error_spotfi","error")
% savefig(fig_error,"Plots/Paper/Localization/triangulation/fig/CDF_error_spotfi")
% print(fig_error,"Plots/Paper/Localization/triangulation/png/CDF_error_spotfi", "-dpng")
% save_PDF_fig(fig_error,"Plots/Paper/Localization/triangulation/png/CDF_error_spotfi")

cd(pwd_str)