close all
clear all

load("mat_files/indoor_map/data_distance_angle_true.mat")
load("mat_files/indoor_map/RX_coordinates.mat")
load("mat_files/indoor_map/points_coordinates.mat")

% Path 1
points_id = [80,21,16,12,8,4,79,28,27,62,70:74,109,110];
% Path 2
points_id = [85, 81, 23, 19, 14, 9, 4, 79, 28, 64, 90, 38, 44, 50, 56];

% Path 4
%points_id = [85,80,21,16,12,8,4,79,28,27,26,68,69,62,63,64,90,91, 92, 93, 41, 46, 51, 56, 61];



load("matlab_scripts/Coordinates/APs_for_Kalman/for_kalman_8.mat")

points_x_traj_real = estimated_point_x_hf_lf_new(points_id);
points_y_traj_real = estimated_point_y_hf_lf_new(points_id);

x = points_x_traj_real(1);
y = points_y_traj_real(1);
initialState = [x;0;y;0];
KF = trackingKF('MotionModel','2D Constant Velocity','State',initialState);

vx = 0.2;
vy = 0.1;
T  = 1.5;
pos = [0:vx*T:2;5:vy*T:6]';

pos = [points_x_traj_real, points_y_traj_real];

figure
hold on
for k = 1:size(points_id,2)
    
    pstates(k,:) = predict(KF, T);
    cstates(k,:) = correct(KF, pos(k,:));
    
end

% Plot
plot(pos(:,1),pos(:,2),'-k.', ...
    cstates(:,1),cstates(:,3),'-o')
xlabel('x [m]')
ylabel('y [m]')
grid
xt  = [x-2 pos(1,1)+0.1 pos(end,1)+0.1];
yt = [y pos(1,2) pos(end,2)];
text(xt,yt,{'First measurement','First position','Last position'})

hold on
plot(points_x(points_id), points_y(points_id), 'gd')
legend('Object position', 'Corrected position')

% CDF
error_corr = sqrt((cstates(:,1)-points_x(points_id)).^2+(cstates(:,3)-points_y(points_id)).^2);
error_jujo = sqrt((points_x_traj_real-points_x(points_id)).^2+(points_y_traj_real-points_y(points_id)).^2);
figure
hold on
cdfplot(error_corr(:))
cdfplot(error_jujo(:))
legend('Corre', 'JUJO')

% Walking speed 1.5 m/s
speed = 1.5;

% We want the error every second
% Path length
length = 0;
time = nan(0, 0);
time(1) = 0;
for id_point=2:size(points_id,2)
    
    prevX = points_x(points_id(id_point-1));
    prevY = points_y(points_id(id_point-1));
    
    currX = points_x(points_id(id_point));
    currY = points_y(points_id(id_point));

    length = length + sqrt((prevX - currX).^2 + (prevY - currY).^2);
    
    time(id_point, 1) = length/speed;
end

figure;
plot(time, error_corr);
hold on
plot(time, error_jujo);

legend("Kalman", "JUJO")
xlabel("Time [s]")
ylabel("Error [m]")

