clear 
close all
clc

load("../../mat_files/indoor/HF/CSI/csi_indoor")
load("../../mat_files/indoor/HF/optimal_rotation_index")
% openfig("../../mat_files/indoor_map/map.fig")
openfig("../../mat_files/indoor_map/map_no_points")

load("../../mat_files/indoor_map/points_coordinates.mat")
load("../../mat_files/indoor_map/RX_coordinates.mat")
load("../../mat_files/indoor_map/data_distance_angle_true.mat")

[n_points, routers, rotations] = size(magnitudes_raw);
power_raw = nan(size(magnitudes_raw));
power_raw_best = nan(n_points,routers);

for ii = 1:n_points
    for jj = 1:routers
        for kk = 1:rotations
           power_ind = magnitudes_raw{ii,jj,kk};
           if (~isempty(power_ind))
               power_raw(ii,jj,kk) = abs(median(median(power_ind,2)))^2;
           end
            
        end
        
        power_raw_best(ii,jj) = power_raw(ii,jj,indoor_optimal_rotation_index(ii,jj));
        
    end
end
    
power_raw_best_db = 10*log10(power_raw_best);

n_colors = 1024;
colormap = jet(n_colors);

power_raw_norm = power_raw_best_db - max(power_raw_best_db);

% do it for router 5


power_router = power_raw_norm(:,5);
grid_power = linspace(0,min(power_router),n_colors);

[~,index_colors] = nanmin(abs(power_router - grid_power),[],2);

points_x_nan = points_x;
points_x_nan(isnan(power_router)) = nan;
points_y_nan = points_y;
points_y_nan(isnan(power_router)) = nan;

scatter(points_x_nan, points_y_nan, [], power_router,'filled')
hold on
scatter(RX_x(5), RX_y(5),'filled', "r")
colorbar



figure, plot(power_router)

power_router_room = power_router([1:25, 75:89]);
figure, plot(power_router_room)

[distance_sort, index_sort] = sort(distances_RX([1:25, 75:89],5));
power_router_sort = power_router_room(index_sort);
figure, plot(power_router_sort)
xticks_variable = xticks;
xticklabels(num2cell(string([0;distance_sort(xticks_variable(2:end))])))
xlabel("Distance to AP [m]")

fc = 60.0e9;
lambda = physconst('LightSpeed')/fc;
% R = 10e3;
L = fspl(distance_sort,lambda)
L_norm = (L - min(L))*(-1);
hold on
plot(L_norm)
ylabel("Attenuation [db]")
legend("Estimated","True")
% legend()
% xticklabels(   