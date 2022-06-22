clear 
clc
close all

pwd_str = pwd;

cd ../../

% create a grid for LOS and NLOS
load("mat_files/indoor_map/points_coordinates.mat")


%% grid for router 10
LOS_RX10 = [true(25,1);false((74-25),1);true(15,1); false((110-89),1)]

openfig("mat_files/indoor_map/map_no_points.fig")
hold on
plot(points_x(LOS_RX10),points_y(LOS_RX10), "*")
plot(points_x(~LOS_RX10),points_y(~LOS_RX10), "*")
title("Router 10")
%% grid for router 4
% LOS_RX4 = [false(1,25),true(1,8),false(1,25),true(1,6),false(1,7),false(1,7)].';
LOS_RX4 = false(110,1);
LOS_RX4([26:31,62:69,105:106]) = true;

openfig("mat_files/indoor_map/map_no_points.fig")
hold on
plot(points_x(LOS_RX4),points_y(LOS_RX4), "*")
plot(points_x(~LOS_RX4),points_y(~LOS_RX4), "*")
title("Router 4")
%% grid for router 5
LOS_RX5 = false(110,1);
LOS_RX5([37:61,90:104 ]) = true;


openfig("mat_files/indoor_map/map_no_points.fig")
hold on
plot(points_x(LOS_RX5),points_y(LOS_RX5), "*")
plot(points_x(~LOS_RX5),points_y(~LOS_RX5), "*")
title("Router 5")

%% grid for router 6
LOS_RX6 = false(110,1);


openfig("mat_files/indoor_map/map_no_points.fig")
hold on
plot(points_x(LOS_RX6),points_y(LOS_RX6), "*")
plot(points_x(~LOS_RX6),points_y(~LOS_RX6), "*")
title("Router 6")

%% grid for router 7
LOS_RX7 = false(110,1);

openfig("mat_files/indoor_map/map_no_points.fig")
hold on
plot(points_x(LOS_RX7),points_y(LOS_RX7), "*")
plot(points_x(~LOS_RX7),points_y(~LOS_RX7), "*")
title("Router 7")

LOS_RX = [LOS_RX4, LOS_RX5, LOS_RX6, LOS_RX7, LOS_RX10];

save("mat_files/indoor_map/Grid_LOS", "LOS_RX")


cd(pwd_str)

