clear
close all
clc

load("mat_files/indoor_map/data_distance_angle_true.mat")
load("mat_files/indoor_map/RX_coordinates.mat")
load("mat_files/indoor_map/points_coordinates.mat")

%fig = openfig("mat_files/indoor_map/map.fig")
% fig2 = openfig("mat_files/indoor_map/map_no_points.fig")
% set(0, 'currentfigure', fig2);
% hold on
% plot(RX_x(2), RX_y(2), "r*")
% plot(RX_x(5), RX_y(5), "r*")

% Path 1
points_id = [80,21,16,12,8,4,79,28,27,62,70:74,109,110];

% Path 2
%points_id = [85, 81, 23, 19, 14, 9, 4, 79, 28, 64, 90, 38, 44, 50, 56];

% Path 3
%points_id = [85,80,21,16,12,8,4,79,28,27,26,68,69,62,63,64,90,91,39, 44, 49, 54, 53, 57, 95, 100];

% Path 4
points_id = [85,80,21,16,12,8,4,79,28,27,26,68,69,62,63,64,90,91, 92, 93, 41, 46, 51, 56, 61];

%points_id = flip([85, 80, 21, 16, 11, 6, 1, 75]);
%points_id = [75, 2, 8, 14, 20]; %Diagonal
%points_id = [27, 64, 90, 38, 44, 50];
%points_id = [85, 81, 23, 19, 15, 10, 5, 79];


% We want the error every second

points_x_traj = points_x(points_id);
points_y_traj = points_y(points_id);

% plot(points_x_traj, points_y_traj, "b*")
% 
% 
% % variance noise of the observation
% power_noise = 0.4;
% 
% traj_x_all_noisy = points_x_traj + randn(size(points_x_traj))*power_noise;
% traj_y_all_noisy = points_y_traj + randn(size(points_y_traj))*power_noise;
% 
% fig3 = openfig("mat_files/indoor_map/map_no_points.fig")
% set(0, 'currentfigure', fig3);
% hold on
% plot(RX_x(2), RX_y(2), "r*")
% plot(RX_x(5), RX_y(5), "r*")
% plot(traj_x_all_noisy, traj_y_all_noisy, "b*")
% 
% 
% % kalman filter
% [traj_x_kalman,traj_y_kalman] = Kalman(traj_x_all_noisy,traj_y_all_noisy);
% fig4 = openfig("mat_files/indoor_map/map_no_points.fig")
% set(0, 'currentfigure', fig4);
% hold on
% plot(RX_x(2), RX_y(2), "r*")
% plot(RX_x(5), RX_y(5), "r*")
% plot(traj_x_kalman, traj_y_kalman, "b*")


load("matlab_scripts/Coordinates/APs_for_Kalman/for_kalman_8.mat")


points_x_traj_real = estimated_point_x_hf_lf_new(points_id);
points_y_traj_real = estimated_point_y_hf_lf_new(points_id);

% fig5 = openfig("mat_files/indoor_map/map_no_points.fig")
% set(0, 'currentfigure', fig5);
% hold on
% plot(RX_x(2), RX_y(2), "r*")
% plot(RX_x(5), RX_y(5), "r*")
% plot(points_x_traj_real, points_y_traj_real, "b*")

% kalman filter
[traj_x_kalman,traj_y_kalman] = Kalman(points_x_traj_real,points_y_traj_real);


XY_kal = SimpleKalman([points_x_traj_real points_y_traj_real]', 1, 1);


fig6 = openfig("mat_files/indoor_map/map_no_points.fig");
set(0, 'currentfigure', fig6);
hold on
plot(RX_x(2), RX_y(2), "r*");
plot(RX_x(5), RX_y(5), "r*");
p1 = plot(points_x_traj, points_y_traj, "-g*");
p2 = plot(points_x_traj_real, points_y_traj_real, "-y*");
p3 = plot(traj_x_kalman, traj_y_kalman, "-b*");
p4 = plot(XY_kal(1, :), XY_kal(2, :), "-ro");

% MSE
hf_error = mean((points_x_traj_real - points_x_traj).^2 + (points_y_traj_real - points_y_traj).^2);
%k1_error = mean((traj_x_kalman - points_x_traj).^2 + (traj_y_kalman - points_y_traj).^2);
k2_error = mean((XY_kal(1, :)' - points_x_traj).^2 + (XY_kal(2, :)' - points_y_traj).^2);

legend([p1,p2,p3, p4],["GT", ['JUJO: ' num2str(hf_error)], "Alex", ['Kalman: ' num2str(k2_error)]]);


% We want the error every second
% Path length
length = 0;
for id_point=2:size(points_id,2)
    
    prevX = points_x(points_id(id_point-1));
    prevY = points_y(points_id(id_point-1));
    
    currX = points_x(points_id(id_point));
    currY = points_y(points_id(id_point));

    length = length + sqrt((prevX - currX).^2 + (prevY - currY).^2);
end

% Walking speed 1.5 m/s
speed = 1.5;

time_it_takes = length/speed;

% We need to interpolate the curves so that they have 1 point per second
[jujo_x_interpolated, jujo_y_interpolated] = interpolateCurve2D(points_x_traj_real, points_y_traj_real, round(time_it_takes));
[kalman_x_interpolated, kalman_y_interpolated] = interpolateCurve2D(XY_kal(1, :)', XY_kal(2, :)', round(time_it_takes));
[gt_x_interpolated, gt_y_interpolated] = interpolateCurve2D(points_x_traj, points_y_traj, round(time_it_takes));

plot(jujo_x_interpolated, jujo_y_interpolated, "-.yx");
plot(kalman_x_interpolated, kalman_y_interpolated, "-.rx");
plot(gt_x_interpolated, gt_y_interpolated, "-.gx");

% Create the graph
hf_error = sqrt((jujo_x_interpolated - gt_x_interpolated).^2 + (jujo_y_interpolated - gt_y_interpolated).^2);
k2_error = sqrt((kalman_x_interpolated - gt_x_interpolated).^2 + (kalman_y_interpolated - gt_y_interpolated).^2);
gt_error = sqrt((gt_x_interpolated - gt_x_interpolated).^2 + (gt_y_interpolated - gt_y_interpolated).^2);

figure
plot([1:size(k2_error, 1)], k2_error);
hold on
plot([1:size(k2_error, 1)], hf_error);
legend("Kalman error", "Jujo error")
xlabel("Time (s)")
ylabel("Error (m)")