clear 
clc
close all



pwd_str = pwd;

cd ../../

addpath("functions")

scenarios = ["indoor", "outdoor", "indoor", "in_out", "imdea"];
scenarios_save = ["indoor", "outdoor", "2ap", "in_out", "imdea"];
rotation_routers_all = {[90, -90, 90, 90, -90],[135, -135, -45, 180],[90, -90, 90, 90, -90],[90, 180],[-90, 90,90]};

routers_all = [5,4,2,2,3];

% for scenario = 1:length(scenarios)
for scenario = 1:length(scenarios)
    error = [];

    load(strcat("mat_files/",scenarios(scenario),"_map/data_distance_angle_true.mat"))
    load(strcat("mat_files/",scenarios(scenario),"_map/points_coordinates.mat"))
    % load(strcat("mat_files/walls.mat")
    load(strcat("mat_files/",scenarios(scenario),"_map/RX_coordinates.mat"))

    load(strcat("mat_files/Jade_all/Data_2D_5_120_Direct_Path_Delay_",scenarios(scenario),".mat"))
%     load("mat_files/Jade/Data_2D_5_120_Direct_Path_Delay_outdoor2.mat")

    routers = routers_all(scenario);
    rotation_routers = rotation_routers_all{scenario};

    % angles_dp = squeeze(aoa_routers(:,:,1));

    if (scenario == 1)
        angles_dp(:,2) = angles_dp(:,2) - 3.5;
        angles_dp(26:38,5) = angles_dp(26:38,5) - 5;

        angles_dp(:,3) = angles_dp(:,3) - 7;
        
        limit_x = 15.3;
        limit_y = 20;

        % move to [x,y] coordinates
        points = [points_x, points_y];
    elseif (scenario == 2)
        angles_dp(:,2) = angles_dp(:,2) + 3;
        angles_dp(:,3) = angles_dp(:,3) - 3;

        limit_x = 17.9;
        limit_y = 24.4;
        points = [points_x.', points_y.'];
        angles_dp(:,[1,4]) = angles_dp(:,[4,1]);

    elseif (scenario == 3)
        angles_dp(:,2) = angles_dp(:,2) - 3.5;
        angles_dp(26:38,5) = angles_dp(26:38,5) - 5;
        angles_dp(:,3) = angles_dp(:,3) - 7;
        
        angles_dp = angles_dp(:,[2,5]);
        routers  = 2;
        rotation_routers = rotation_routers(:,[2,5]);
        RX_x = RX_x(:,[2,5]);
        RX_y = RX_y(:,[2,5]);
        limit_x = 15.3;
        limit_y = 20;

        points = [points_x, points_y];
    elseif (scenario == 4)

        limit_x = 25;
        limit_y = 18;
        points_x = points_x + 7;
        points_y = points_y + 2;
        points = [points_x, points_y];
        RX_x = RX_x + 7;
        RX_y = RX_y + 2;
%         angles_dp(:,[1 2]) = angles_dp(:,[2 1])

%         angles_dp = angles_RX;
    elseif(scenario == 5)
        limit_x = 12;
        limit_y = 25;
        angles_dp(:,3) = angles_dp(:,3) - 6;
        angles_dp(:,2) = angles_dp(:,2) + 6;
        angles_dp([19,20,27:31],1) = angles_dp([19,20,27:31],1)*(-1);
        points = [points_x, points_y];
    end
    

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




    % put the max and min number to create the grid
    min_angle = -90;
    max_angle = +90;

    % inverse the symbol
    % angles_dp = (-1)*angles_dp;

    % Plot_Map_No_Points()



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
        heatmap_total = nan([size(X),routers]);

        % loop to iterate over aps
        for id_router = 1:routers


            % take the angles to test
            angles_test = min_angle:0.1:max_angle;

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

    %         % take the first decimal
    %         angles_ap = round(rad2deg(atan(diff_X./diff_Y)),1); 
             % take the first decimal
            angles_ap = round(rad2deg(angle(diff_X + 1i*diff_Y)),1);    
            angles_ap = angles_ap - rotation_routers(id_router);
            angles_ap(angles_ap < -180) = angles_ap(angles_ap < -180) + 360;
            angles_ap = round(angles_ap,1);

            angles_test = round(angles_test,1);


            % loop to assing to each angle to searcht its power
            for id_angle = 1:length(angles_test)

                index_change = angles_ap == angles_test(id_angle);
                heatmap(index_change) = power_dist(id_angle);

            end

            % if the heatmap is 0 means that there is no power, but it is
            % better to set as the th
            index_0 = heatmap == 0;
            heatmap(index_0) = th;
% 
%             figure, pcolor(heatmap)
%             colorbar

            % save to a big matriz
            heatmap_total(:,:,id_router) = heatmap;
        end


        % take the product of the heatmaps
        heatmap = prod(heatmap_total,3);
% 
%         figure, pcolor(heatmap)
%         colorbar

        % take the maximum
        [max_num,max_idx] = max(heatmap(:));  
        % find it in the matrix
        [row, column] = find(heatmap == max_num);
        % take the position
        solution = [X(row(1), column(1)),Y(row(1), column(1))];
        solution_total(id_point,:) = solution;

        % estimate the error
        error(id_point,1) = sqrt(sum(((solution) - position).^2,2));
    end
    [~, index_max] = max(error);
    % error(index_max) = 10;
    fig_error = figure;
    cdfplot(error);
    save(strcat("mat_files/CDF_data/",scenarios_save(scenario),"/error_spotfi"),"error")
    % savefig(fig_error,"Plots/Paper/Localization/triangulation/fig/CDF_error_spotfi")
    % print(fig_error,"Plots/Paper/Localization/triangulation/png/CDF_error_spotfi", "-dpng")
    % save_PDF_fig(fig_error,"Plots/Paper/Localization/triangulation/png/CDF_error_spotfi")
end
cd(pwd_str)