close all
clear all
clc

%%

% Load the ground truth
load('mat_files/indoor_map/RX_coordinates.mat');
load('mat_files/indoor_map/points_coordinates.mat');
load('mat_files/indoor_map/data_distance_angle_true.mat');

% LF data
load("mat_files/mD_track/LF/Data_3D.mat")
load('mat_files/FTM/LF/estimated_distance.mat')

% HF data
load('mat_files/indoor/HF/CSI/csi_indoor.mat')
load('mat_files/indoor/HF/FTM/ftm_indoor.mat')

% Coordinates
load('mat_files/Coordinates/data_coordinates_hf.mat')
load('mat_files/Coordinates/data_coordinates_lf.mat')

% Offsets
load("mat_files/indoor/HF/CSI/aoa_offset.mat")

% We want to localize based on the multipath components and ToF
%% Estimate the positions
estimated_points_x_hf = cell(110, 5, 8);
estimated_points_y_hf = cell(110, 5, 8);

for id_point = 1:110
    for id_router = 1:5
        for id_rotation = 1:8

            % For each path
            paths = azimut_raw{id_point, id_router, id_rotation};
            
            path_x_aux = nan(0, 0);
            path_y_aux = nan(0, 0);
            
            % If there are no paths, continue
            if size(paths, 2) == 0
               continue 
            end
            
            for id_path=1:size(paths, 2)
            
                estimated_point_x_aux = sind((paths(id_path) - aoa_offset(id_router))*-1)*calculated_distance(id_point,id_router, id_rotation);
                estimated_point_y_aux = cosd((paths(id_path) - aoa_offset(id_router))*-1)*calculated_distance(id_point,id_router, id_rotation);

        %         figure, plot(estimated_point_x_p(id_point,id_router), estimated_point_y_p(id_point,id_router), "*")

                if (id_router == 1 || id_router == 3 || id_router == 4)
                    path_x_aux(end+1) = RX_x(id_router) - (estimated_point_x_aux);
                    path_y_aux(end+1) = RX_y(id_router) - (estimated_point_y_aux*(-1));
                else
                    path_x_aux(end+1) = RX_x(id_router) - (estimated_point_x_aux*(-1));
                    path_y_aux(end+1) = RX_y(id_router) - estimated_point_y_aux;
                end
            end
            
            estimated_points_x_hf{id_point, id_router, id_rotation} = path_x_aux;
            estimated_points_y_hf{id_point, id_router, id_rotation} = path_y_aux;
        end
    end
end

save("mat_files/Angle_profile/data_coordinates_hf_all_paths","estimated_points_x_hf", "estimated_points_y_hf");

% LF
estimated_points_x_lf = cell(110, 5);
estimated_points_y_lf = cell(110, 5);

for id_point = 1:110
    for id_router = 1:5

        paths = AoA{id_point, id_router};
        paths = paths{1};
        
        path_x_aux = nan(0, 0);
        path_y_aux = nan(0, 0);
            
        % Remove the offsets
        if id_point > 1 && id_point < 75 && id_router == 2
            
            paths = paths - 3.5;
        end
        
        if id_point > 1 && id_point < 75 && id_router == 3
            
            paths = paths - 7;
        end
        
        if id_point > 25 && id_point < 39 && id_router == 5
            
            paths = paths - 5;
        end
        
        % If there are no paths, continue
        if size(paths, 2) == 0
           continue 
        end
        
        for id_path=1:size(paths, 2)

            estimated_point_x_p(id_point,id_router) = sind(paths(id_path))*(estimated_distance(id_point,id_router));
            estimated_point_y_p(id_point,id_router) = cosd(paths(id_path))*(estimated_distance(id_point,id_router));

            if (id_router == 1 || id_router == 3 || id_router == 4)
                path_x_aux(end+1) = RX_x(id_router) - (estimated_point_x_p(id_point,id_router));
                path_y_aux(end+1) = RX_y(id_router) - (estimated_point_y_p(id_point,id_router)*(-1));
            else
                path_x_aux(end+1) = RX_x(id_router) - (estimated_point_x_p(id_point,id_router)*(-1));
                path_y_aux(end+1) = RX_y(id_router) - estimated_point_y_p(id_point,id_router);
            end
        end
        
        estimated_points_x_lf{id_point, id_router} = path_x_aux;
        estimated_points_y_lf{id_point, id_router} = path_y_aux;
    end
end

save("mat_files/Angle_profile/data_coordinates_lf_all_paths","estimated_points_x_lf", "estimated_points_y_lf");

%% Draw them
% Draw only the optimal rotation for HF
%load('mat_files/indoor/HF/optimal_rotation_index.mat');

% Chose the antenna array from LF
load('mat_files/Coordinates/chosen_rotation_with_lf.mat')


chosen_points = [28];
chosen_routers = [1:5];

% Generate a different color for each point
N = size(chosen_routers, 2);
fun = @(m)sRGB_to_OSAUCS(m,true,true); % recommended OSA-UCS
rgb = maxdistcolor(N,fun);

% Get a different marker for each AP
markers = {'d','o','*','s','x','.','+','^','v','>','<','p','h'};

for rotation_id=1:4
    
    uiopen('mat_files/indoor_map/map.fig',1);
    hold on
    pause(1)

    plot(points_x, points_y, "*", "LineWidth", 4, 'color', [1, 1, 1]);
    plot(RX_x, RX_y, "*", "LineWidth", 4, 'color', [1, 1, 1]);
    
    for point_id=1:size(chosen_points, 2)
        for router_id=1:size(chosen_routers, 2)
        
            scatter(RX_x(router_id), RX_y(router_id), 50, "*", "LineWidth", 4,'MarkerEdgeColor', rgb(router_id, :), 'MarkerFaceColor', rgb(router_id, :));

            chosen_point = chosen_points(point_id);
            chosen_router = chosen_routers(router_id);
            
            % Chose array 1 or 2
            if chosen_rotation_with_lf(chosen_point, chosen_router, rotation_id) == 1
                chosen_rotation = rotation_id;
            else
                chosen_rotation = rotation_id+4;
            end
            
            hold on

            if ~isnan(chosen_rotation_with_lf(chosen_point, chosen_router))
                x = estimated_points_x_hf{chosen_point, chosen_router, chosen_rotation};
                y = estimated_points_y_hf{chosen_point, chosen_router, chosen_rotation};

                scatter(x, y, 100, 'o', 'MarkerEdgeColor', rgb(router_id, :));

                % Draw a circle around the chosen AP
                circle(RX_x(chosen_router), RX_y(chosen_router), calculated_distance(chosen_point, chosen_router, chosen_rotation), '-', rgb(router_id, :)) %HF

            end

            x = estimated_points_x_lf{chosen_point, chosen_router};
            y = estimated_points_y_lf{chosen_point, chosen_router};

            scatter(x, y, 100, 'x', 'MarkerEdgeColor', rgb(router_id, :));

            % The estimation point
            plot(points_x(chosen_point), points_y(chosen_point), "*", "LineWidth", 4, 'color', 'black');

            % Draw a circle around the chosen AP
            circle(RX_x(chosen_router), RX_y(chosen_router), estimated_distance(chosen_point, chosen_router), '-.', rgb(router_id, :)) %LF

        end
    end
    
    hold off
    pbaspect([1 1 1]);
    
    title(['Rotation ' num2str(rotation_id)])
end

%% Everything at the same time
% Draw only the optimal rotation for HF
%load('mat_files/indoor/HF/optimal_rotation_index.mat');

% Chose the antenna array from LF
load('mat_files/Coordinates/chosen_rotation_with_lf.mat')

chosen_points = [28];
chosen_routers = [1:5];

% Generate a different color for each point
N = size(chosen_routers, 2);
fun = @(m)sRGB_to_OSAUCS(m,true,true); % recommended OSA-UCS
rgb = maxdistcolor(N,fun);

% Get a different marker for each AP
markers = {'d','o','*','s','x','.','+','^','v','>','<','p','h'};

uiopen('mat_files/indoor_map/map.fig',1);
hold on
pause(1)

plot(points_x, points_y, "*", "LineWidth", 4, 'color', [1, 1, 1]);
plot(RX_x, RX_y, "*", "LineWidth", 4, 'color', [1, 1, 1]);

    
for point_id=1:size(chosen_points, 2)
    for router_id=1:size(chosen_routers, 2)
        
        scatter(RX_x(router_id), RX_y(router_id), 50, "*", "LineWidth", 4,'MarkerEdgeColor', rgb(router_id, :), 'MarkerFaceColor', rgb(router_id, :));

        chosen_point = chosen_points(point_id);
        chosen_router = chosen_routers(router_id);
            
        x = estimated_points_x_lf{chosen_point, chosen_router};
        y = estimated_points_y_lf{chosen_point, chosen_router};

        scatter(x, y, 100, 'x', 'MarkerEdgeColor', rgb(router_id, :));

        % The estimation point
        plot(points_x(chosen_point), points_y(chosen_point), "*", "LineWidth", 4, 'color', 'black');

        % Draw a circle around the chosen AP
        circle(RX_x(chosen_router), RX_y(chosen_router), estimated_distance(chosen_point, chosen_router), '-.', rgb(router_id, :)) %LF
            
        for rotation_id=1:4

            % Chose array 1 or 2
            if chosen_rotation_with_lf(chosen_point, chosen_router, rotation_id) == 1
                chosen_rotation = rotation_id;
            else
                chosen_rotation = rotation_id+4;
            end
            
            hold on

            if ~isnan(chosen_rotation_with_lf(chosen_point, chosen_router))
                x = estimated_points_x_hf{chosen_point, chosen_router, chosen_rotation};
                y = estimated_points_y_hf{chosen_point, chosen_router, chosen_rotation};

                scatter(x, y, 100, markers{rotation_id}, 'MarkerEdgeColor', rgb(router_id, :));

                % Draw a circle around the chosen AP
                circle(RX_x(chosen_router), RX_y(chosen_router), calculated_distance(chosen_point, chosen_router, chosen_rotation), '-', rgb(router_id, :)) %HF

            end
        end
    end 
end

hold off
pbaspect([1 1 1]);

title(['All rotations'])
