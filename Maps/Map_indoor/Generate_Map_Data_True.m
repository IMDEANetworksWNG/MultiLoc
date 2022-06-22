clear
close all
clc

load("walls")
load("RX_coordinates.mat")
load("points_coordinates.mat")
load("data_distance_angle_true.mat")

mkdir("../../mat_files")
mkdir("../../mat_files/indoor_map")

figure,
hold on
% prx_51 = plot(RX_x_1W51, RX_y_1W51, "r*", "LineWidth", 4)
% prx_47_ = plot(RX_x_1W47, RX_y_1W47, "r*", "LineWidth", 4)
% prx_54 = plot(RX_x_1W54, RX_y_1W54, "r*", "LineWidth", 4)
% prx_58 = plot(RX_x_1W58, RX_y_1W58, "r*", "LineWidth", 4)
% 
% prx_corridor = plot(RX_x_corridor, RX_y_corridor, "r*", "LineWidth", 4)

p_wall_1W51 = plot(walls_rigth_x_1W51, walls_rigth_y_1W51, "k", "LineWidth", 1)
plot(walls_left_x_1W51, walls_left_y_1W51, "k", "LineWidth", 1)
plot(walls_bottom_x_1W51, walls_bottom_y_1W51, "k", "LineWidth", 1)
plot(walls_top_x_1W51, walls_top_y_1W51, "k", "LineWidth", 1)

p_wall_1W49 = plot(walls_rigth_x_1W49, walls_rigth_y_1W49, "k", "LineWidth", 1)
plot(walls_left_x_1W49, walls_left_y_1W49, "k", "LineWidth", 1)
plot(walls_bottom_x_1W49, walls_bottom_y_1W49, "k", "LineWidth", 1)
plot(walls_top_x_1W49, walls_top_y_1W49, "k", "LineWidth", 1)

p_wall_1W47 = plot(walls_rigth_x_1W47, walls_rigth_y_1W47, "k", "LineWidth", 1)
plot(walls_left_x_1W47, walls_left_y_1W47, "k", "LineWidth", 1)
plot(walls_bottom_x_1W47, walls_bottom_y_1W47, "k", "LineWidth", 1)
plot(walls_top_x_1W47, walls_top_y_1W47, "k", "LineWidth", 1)

p_wall_1W58 = plot(walls_rigth_x_1W58, walls_rigth_y_1W58, "k", "LineWidth", 1)
hold on
plot(walls_left_x_1W58, walls_left_y_1W58, "k", "LineWidth", 1)
plot(walls_bottom_x_1W58, walls_bottom_y_1W58, "k", "LineWidth", 1)
plot(walls_top_x_1W58, walls_top_y_1W58, "k", "LineWidth", 1)

p_wall_1W56 = plot(walls_rigth_x_1W56, walls_rigth_y_1W56, "k", "LineWidth", 1)
plot(walls_left_x_1W56, walls_left_y_1W56, "k", "LineWidth", 1)
plot(walls_bottom_x_1W56, walls_bottom_y_1W56, "k", "LineWidth", 1)
plot(walls_top_x_1W56, walls_top_y_1W56, "k", "LineWidth", 1)

p_wall_1W54 = plot(walls_rigth_x_1W54, walls_rigth_y_1W54, "k", "LineWidth", 1)
plot(walls_left_x_1W54, walls_left_y_1W54, "k", "LineWidth", 1)
plot(walls_bottom_x_1W54, walls_bottom_y_1W54, "k", "LineWidth", 1)
plot(walls_top_x_1W54, walls_top_y_1W54, "k", "LineWidth", 1)

savefig("../../mat_files/indoor_map/map_no_points")


p_points = plot(points_x,points_y, "b*", "LineWidth", 4);
p_routers = plot(RX_x,RX_y, "r*", "LineWidth", 4);


AP_labels_HF = [45,44,47,46,43];
AP_labels_LF = [4,5,6,7,10];

text(RX_x + 0.1, RX_y + 0.1, strcat(string(AP_labels_LF), "/", string(AP_labels_HF)));

p_points = plot(points_x,points_y, "b*", "LineWidth", 4);
text(points_x + 0.1, points_y + 0.1, string(1:110));

ylim([-0.5 (walls_top_y_1W58(end) + 0.5)])
xlim([-0.5 (walls_rigth_x_1W58(1) + 0.5)])

savefig("../../mat_files/indoor_map/map")


% take the distances between APs to points
distances_RX = sqrt((points_x - RX_x).^2 + (points_y - RX_y).^2);

% take the angles between APs to points
diff_RX_x = points_x - RX_x;
diff_RX_y = points_y - RX_y;
angles_RX = rad2deg(atan(diff_RX_x./diff_RX_y));

% Point labels
point_labels = 1:110;



save("../../mat_files/indoor_map/points_coordinates", "points_x", "points_y", "point_labels")
save("../../mat_files/indoor_map/RX_coordinates", "RX_x","RX_y")
save("../../mat_files/indoor_map/data_distance_angle_true", "angles_RX", "distances_RX", "AP_labels_LF", "AP_labels_HF")

clearvars -except walls_*
% clearvars walls_x*
% clearvars walls_y*

save("../../mat_files/indoor_map/walls")
