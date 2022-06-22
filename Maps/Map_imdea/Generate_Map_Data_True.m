close all
clc
clear

mkdir ("../../mat_files/")
mkdir ("../../mat_files/imdea_map")

%% Open area 2E


offset_x_izq_2E_open = 0.36;
offset_x_drch_2E_open = 0.36;

limit_x_2E_open = offset_x_izq_2E_open + 11 * 0.6 + offset_x_drch_2E_open;
limit_y_2E_open = 30 * 0.6;

RX_x_2E_open_1 = [1]*0.6 + offset_x_izq_2E_open;
RX_y_2E_open_1 = [1]*0.6;

RX_x_2E_open_2 = [9]*0.6 + offset_x_izq_2E_open;
RX_y_2E_open_2 = [24]*0.6;

points_y_2E_open = repmat(5:3:21,2,1);
points_x_2E_open = ones(size(points_y_2E_open)).*[1,10].';

points_y_2E_open = points_y_2E_open(:);
points_x_2E_open = points_x_2E_open(:);

points_y_2E_open_desk = repmat([9 12 17],2,1);
points_x_2E_open_desk = ones(size(points_y_2E_open_desk)).*[3 8].';

points_y_2E_open_desk = points_y_2E_open_desk(:);
points_x_2E_open_desk = points_x_2E_open_desk(:);

points_y_2E_office = [23 + 4;23 + 4];
points_x_2E_office = [1; -2];

points_y_2E_office2 = [23 + 4;23 + 4];
points_x_2E_office2 = [11; 14];

% points_y_2E_office = points_y_2E_office * 0.6;
% points_x_2E_office = (points_x_2E_office * 0.6) + offset_x_izq_2E_open;

points_y_2E_open = [points_y_2E_open;points_y_2E_open_desk;points_y_2E_office];
points_x_2E_open = [points_x_2E_open;points_x_2E_open_desk;points_x_2E_office];


points_y_2E_open = points_y_2E_open * 0.6;
points_x_2E_open = (points_x_2E_open * 0.6) + offset_x_izq_2E_open;

points_y_2E_office2 = points_y_2E_office2 * 0.6;
points_x_2E_office2 = (points_x_2E_office2 * 0.6) + offset_x_izq_2E_open;

figure('units','normalized','outerposition',[0 0 0.4 0.8])
prx1 = plot(RX_x_2E_open_1,RX_y_2E_open_1, "*", "LineWidth", 4)
hold on
prx2 = plot(RX_x_2E_open_2, RX_y_2E_open_2, "*", "LineWidth", 4)
p_points = plot(points_x_2E_open,points_y_2E_open, "*", "LineWidth", 4)


% corridor
% space between end of open arean and bathroom
% bathroom and space between bathroom and corridor
% limit_y_2E_corridor = 0.6 + (4.08 + 0.115) + 0.38 - 0.2;
limit_y_2E_corridor = (9*0.5) + 0.15 + (1.05 - 0.6);
limit_x_2E_corridor = 0;

limit_x_2E_open = limit_x_2E_open + limit_x_2E_corridor;
limit_y_2E_open = limit_y_2E_open + limit_y_2E_corridor;

points_y_2E_corridor = (30 * 0.6 + (1.05 - 0.6));
points_y_2E_corridor = points_y_2E_corridor:1.5:(points_y_2E_corridor + 1.5*2)
points_x_2E_corridor = (offset_x_izq_2E_open + 11 * 0.6 + offset_x_drch_2E_open) - (4 * 0.6 + 0.48);
points_x_2E_corridor = ones(1, length(points_y_2E_corridor)) * points_x_2E_corridor;
p_points = plot(points_x_2E_corridor,points_y_2E_corridor, "*", "LineWidth", 4)

xlim([-0.5 (limit_x_2E_open + 0.5)]);
ylim([-0.5 (limit_y_2E_open + 0.5)]);

%% Marconi

offset_x_izq_2E_marconi = 0.2;
offset_x_drch_2E_marconi = 0.488;

wall_x_offset_marconi = 0.6 - (offset_x_izq_2E_marconi + offset_x_drch_2E_open);

limit_x_2E_marconi = wall_x_offset_marconi + offset_x_izq_2E_marconi + (6 * 0.6) + offset_x_drch_2E_marconi;
limit_y_2E_marconi = 12 * 0.6;

points_y_marconi = repmat(5:3:11,2,1);
points_x_marconi = ones(size(points_y_marconi)).*[1,5].';

points_y_marconi = points_y_marconi(:);
points_x_marconi = points_x_marconi(:);

points_y_marconi = points_y_marconi * 0.6;
points_x_marconi = wall_x_offset_marconi + offset_x_izq_2E_marconi + (points_x_marconi * 0.6) + offset_x_drch_2E_marconi;

rx_x_marconi = 3;
rx_y_marconi = 1;

rx_y_marconi = rx_y_marconi * 0.6;
rx_x_marconi = wall_x_offset_marconi + offset_x_izq_2E_marconi + (rx_x_marconi * 0.6) + offset_x_drch_2E_marconi;


figure('units','normalized','outerposition',[0 0 0.4 0.8])
p_points = plot(points_x_marconi,points_y_marconi, "*", "LineWidth", 4)
hold on
p_points = plot(rx_x_marconi,rx_y_marconi, "*", "LineWidth", 4)

% hold on
% prx = plot(RX_x_2E_open_2, RX_y_2E_open_2, "*", "LineWidth", 4)


%% Join open area 2S with Marconi

points_x_marconi = points_x_marconi + limit_x_2E_open;
points_y_marconi = points_y_marconi + 0;

rx_x_marconi = rx_x_marconi + limit_x_2E_open;
rx_y_marconi = rx_y_marconi + 0;

points_x = [points_x_2E_open;points_x_marconi;points_x_2E_corridor.';points_x_2E_office2];
points_y = [points_y_2E_open;points_y_marconi;points_y_2E_corridor.';points_y_2E_office2];

% [43 52 54]
RX_x = [RX_x_2E_open_2;rx_x_marconi;RX_x_2E_open_1];
RX_y = [RX_y_2E_open_2;rx_y_marconi;RX_y_2E_open_1];

figure('units','normalized','outerposition',[0 0 0.4 0.8])
% p_points = plot(RX_x_2E_open_1,RX_y_2E_open_1, "*", "LineWidth", 4)
hold on
% prx = plot(RX_x_2E_open_2, RX_y_2E_open_2, "*", "LineWidth", 4)
% prx = plot(rx_x_marconi, rx_y_marconi, "*", "LineWidth", 4)
prx = plot(RX_x, RX_y, "*", "LineWidth", 4)


% p_points = plot(points_x_2E_open,points_y_2E_open, "*", "LineWidth", 4)
% p_points = plot(points_x_marconi,points_y_marconi, "*", "LineWidth", 4)
p_points = plot(points_x,points_y, "*", "LineWidth", 4)

text(points_x + 0.2, points_y + 0.2,string(1:length(points_x)))
% aps = miktrotik, asus CSI and asus FTM
aps_string = ["43-222-2", "52-223-4", "54-10-1"];
text(RX_x + 0.2, RX_y + 0.2,aps_string)

distances_RX = sqrt((points_x - RX_x.').^2 + (points_y - RX_y.').^2);
diff_RX_x = points_x - RX_x.';
diff_RX_y = points_y - RX_y.';
angles_RX = rad2deg(atan(diff_RX_x./diff_RX_y));


save("../../mat_files/imdea_map/data_distance_angle_true", "distances_RX", "angles_RX")
save("../../mat_files/imdea_map/points_coordinates", "points_x", "points_y")
save("../../mat_files/imdea_map/RX_coordinates", "RX_x", "RX_y")

savefig(figure(3), "../../mat_files/imdea_map/map")
% xlim([-0.5 (limit_x_2E_open + limit_x_2E_marconi + 0.5)]);
% ylim([-0.5 (limit_y_2E_open + 0.5)]);

% %% Corridor
% 
% offset_x_izq_corridor = 0.23;
% offset_x_drch_corridor = 0.35;
% 
% offset_y_below_corridor = 0.41;
% offset_y_above_corridor = 0.31 + 1.57;
% 
% 
% 
% 
% limit_x_corridor = offset_x_izq_corridor + 11 * 0.5 + offset_x_drch_corridor;
% limit_y_corridor = offset_y_below_corridor + 34 * 0.5 + offset_y_above_corridor;
% 
% points_y_corridor = 4:3:31;
% points_x_corridor = ones(1,length(points_y_corridor))*10;
% 
% points_x_corridor = offset_x_izq_corridor + points_x_corridor * 0.5;
% points_y_corridor = offset_y_below_corridor + points_y_corridor * 0.5;
% 
% RX_y_corridor = 4;
% RX_x_corridor = 2;
% 
% RX_x_corridor = offset_x_izq_corridor + RX_x_corridor * 0.5;
% RX_y_corridor = offset_y_below_corridor + RX_y_corridor * 0.5;
% 
% figure('units','normalized','outerposition',[0 0 0.4 0.8])
% p_points = plot(RX_x_corridor,RX_y_corridor, "*", "LineWidth", 4)
% hold on
% p_points = plot(points_x_corridor,points_y_corridor, "*", "LineWidth", 4)
% 
% %% Join corridor with open area 2E
% distance_y_corridor_2E_open = limit_y_corridor - limit_x_2E_open + offset_x_izq_2E_open;
% distance_x_corridor_2E_open = limit_x_corridor + limit_y_2E_open;
% 
% RX_y_2E_open_1_aux = RX_x_2E_open_1 + distance_y_corridor_2E_open;
% RX_x_2E_open_1 = -RX_y_2E_open_1 + distance_x_corridor_2E_open;
% RX_y_2E_open_1 = RX_y_2E_open_1_aux;
% 
% RX_y_2E_open_2_aux = RX_x_2E_open_2 + distance_y_corridor_2E_open;
% RX_x_2E_open_2 = -RX_y_2E_open_2 + distance_x_corridor_2E_open;
% RX_y_2E_open_2 = RX_y_2E_open_2_aux;
% 
% points_y_2E_open_aux = points_x_2E_open + distance_y_corridor_2E_open;
% points_x_2E_open = -points_y_2E_open + distance_x_corridor_2E_open;
% points_y_2E_open = points_y_2E_open_aux;
% 
% points_y_2E_office_aux = points_x_2E_office + distance_y_corridor_2E_open;
% points_x_2E_office = -points_y_2E_office + distance_x_corridor_2E_open;
% points_y_2E_office = points_y_2E_office_aux;
% 
% points_y_2E_corridor_aux = points_x_2E_corridor + distance_y_corridor_2E_open;
% points_x_2E_corridor = -points_y_2E_corridor + distance_x_corridor_2E_open;
% points_y_2E_corridor = points_y_2E_corridor_aux;
% 
% 
% 
% points_y_marconi_aux = points_x_marconi + distance_y_corridor_2E_open;
% points_x_marconi = -points_y_marconi + distance_x_corridor_2E_open;
% points_y_marconi = points_y_marconi_aux;
% 
% figure('units','normalized','outerposition',[0 0 0.4 0.8])
% prx1 = plot(RX_x_2E_open_1,RX_y_2E_open_1, "*", "LineWidth", 4)
% hold on
% prx2 = plot(RX_x_2E_open_2, RX_y_2E_open_2, "*", "LineWidth", 4)
% p_points = plot(RX_x_corridor,RX_y_corridor, "*", "LineWidth", 4)
% p_points = plot(points_x_2E_open,points_y_2E_open, "*", "LineWidth", 4)
% p_points = plot(points_x_corridor,points_y_corridor, "*", "LineWidth", 4)
% p_points = plot(points_x_marconi,points_y_marconi, "*", "LineWidth", 4)
% p_points = plot(points_x_2E_corridor,points_y_2E_corridor, "*", "LineWidth", 4)
% 
% 
% %% Open area 2S
% 
% offset_x_izq_2S_open = 0.12;
% offset_x_drch_2S_open = 0.12;
% 
% limit_x_2S_open = offset_x_izq_2S_open + 12 * 0.6 + offset_x_drch_2S_open;
% limit_y_2S_open = 30 * 0.6;
% 
% points_y_2S_open = repmat(5:3:23,2,1);
% points_x_2S_open = ones(size(points_y_2S_open)).*[1,11].';
% 
% points_y_2S_offices1 = [27,27];
% points_x_2S_offices1 = [13 16];
% 
% points_y_2S_offices2 = [22,22];
% points_x_2S_offices2 = 0 - [1,4];
% 
% points_y_2S_open = points_y_2S_open(:);
% points_x_2S_open = points_x_2S_open(:);
% 
% points_y_2S_offices2 = points_y_2S_offices2(:);
% points_x_2S_offices2 = points_x_2S_offices2(:);
% 
% points_y_2S_offices1 = points_y_2S_offices1(:);
% points_x_2S_offices1 = points_x_2S_offices1(:);
% 
% points_y_2S_open_desk = repmat([9 13 17],2,1);
% points_x_2S_open_desk = ones(size(points_y_2S_open_desk)).*[3 8].';
% 
% points_y_2S_open_desk = points_y_2S_open_desk(:);
% points_x_2S_open_desk = points_x_2S_open_desk(:);
% 
% 
% points_y_2S_open = [points_y_2S_open;points_y_2S_offices1;points_y_2S_offices2;points_y_2S_open_desk];
% points_x_2S_open = [points_x_2S_open;points_x_2S_offices1;points_x_2S_offices2;points_x_2S_open_desk];
% 
% points_y_2S_open = points_y_2S_open * 0.6;
% points_x_2S_open = points_x_2S_open * 0.6 + offset_x_izq_2S_open;
% 
% % points_y_2S_open = [points_y_2S_open;*0.6];
% % points_x_2S_open = [points_x_2S_open;points_x_2S_offices2*0.6];
% 
% RX_x_2S_open_1 = [6]*0.6 + offset_x_izq_2S_open;
% RX_y_2S_open_1 = [1]*0.6;
% 
% RX_x_2S_open_2 = [10]*0.6 + offset_x_izq_2S_open;
% RX_y_2S_open_2 = [27]*0.6;
% 
% 
% figure('units','normalized','outerposition',[0 0 0.4 0.8])
% p_points = plot(RX_x_2S_open_1,RX_y_2S_open_1, "r*", "LineWidth", 4)
% hold on
% prx = plot(RX_x_2S_open_2, RX_y_2S_open_2, "r*", "LineWidth", 4)
% p_points = plot(points_x_2S_open, points_y_2S_open, "b*", "LineWidth", 4)
% 
% % xlim([-0.5 (limit_x_2S_open + 0.5)]);
% % ylim([-0.5 (limit_y_2S_open + 0.5)]);
% 
% 
% 
% 
% 
% %% Join open area everything
% 
% total_area_x = 6.113 + 11 * 0.4 + 0.15 + 4.653;
% diff_corridor = 0.59;
% diff_area = 0.47;
% 
% limit_x_corridor - diff_corridor
% end_corridor_area = 4.653 + 0.47 - 0.59;
% start_corridor_area_x = total_area_x - limit_x_corridor - end_corridor_area;
% 
% offices_2S = (total_area_x - limit_x_2S_open)/2;
% 
% 
% 
% distance_x_2S_corridor = start_corridor_area_x - (offices_2S + offset_x_izq_2S_open)
% distance_y_2S_corridor = 0.1*2 + 0.17 + 0.17 + 11*0.4 + limit_y_2S_open
% 
% %% area
% distance_x_2S_area = total_area_x - offices_2S -offset_x_izq_2S_open - 4.653;
% distance_y_2S_area = limit_y_2S_open + 0.1;
% 
% point_x_area = [2 6 6 10];
% point_y_area = [7 7 2 2];
% 
% point_x_area = point_x_area * 0.4 + 0.15;
% point_y_area = point_y_area * 0.4 + 0.17;
% 
% point_y_area_aux = point_x_area + distance_y_2S_area;
% point_x_area = -point_y_area + distance_x_2S_area
% point_y_area = point_y_area_aux;
% 
% 
% 
% % 2E
% RX_x_2E_open_1 = RX_x_2E_open_1 + distance_x_2S_corridor;
% RX_y_2E_open_1 = RX_y_2E_open_1 + distance_y_2S_corridor;
% 
% RX_x_2E_open_2 = RX_x_2E_open_2 + distance_x_2S_corridor;
% RX_y_2E_open_2 = RX_y_2E_open_2 + distance_y_2S_corridor;
% 
% RX_x_corridor = RX_x_corridor + distance_x_2S_corridor;
% RX_y_corridor = RX_y_corridor + distance_y_2S_corridor;
% 
% 
% points_x_2E_open = points_x_2E_open + distance_x_2S_corridor;
% points_y_2E_open = points_y_2E_open + distance_y_2S_corridor;
% 
% points_x_2E_office = points_x_2E_office + distance_x_2S_corridor;
% points_y_2E_office = points_y_2E_office + distance_y_2S_corridor;
% 
% 
% points_x_corridor = points_x_corridor + distance_x_2S_corridor;
% points_y_corridor = points_y_corridor + distance_y_2S_corridor;
% 
% points_x_marconi = points_x_marconi + distance_x_2S_corridor;
% points_y_marconi = points_y_marconi + distance_y_2S_corridor;
% 
% points_x_2E_corridor = points_x_2E_corridor + distance_x_2S_corridor;
% points_y_2E_corridor = points_y_2E_corridor + distance_y_2S_corridor;
% 
% prx1 = plot(RX_x_2E_open_1,RX_y_2E_open_1, "r*", "LineWidth", 4)
% prx2 = plot(RX_x_2E_open_2, RX_y_2E_open_2, "r*", "LineWidth", 4)
% p_points = plot(RX_x_corridor,RX_y_corridor, "r*", "LineWidth", 4)
% p_points = plot(points_x_2E_open,points_y_2E_open, "b*", "LineWidth", 4)
% p_points = plot(points_x_corridor,points_y_corridor, "b*", "LineWidth", 4)
% p_points = plot(points_x_marconi,points_y_marconi, "b*", "LineWidth", 4)
% p_points = plot(points_x_2E_corridor,points_y_2E_corridor, "b*", "LineWidth", 4)
% p_points = plot(point_x_area,point_y_area, "b*", "LineWidth", 4)
% p_points = plot(points_x_2E_office,points_y_2E_office, "b*", "LineWidth", 4)
% 
% 
% 
% limit_y = distance_y_2S_corridor + distance_y_corridor_2E_open + limit_x_2E_marconi + limit_x_2E_open
% limit_x = distance_x_2S_corridor + distance_x_corridor_2E_open
% 
% % xlim([-5 (limit_x + 0.5)])
% % ylim([-0.5 (limit_y + 0.5)])
% 
% fig = figure(6);
% fig.CurrentAxes
% ax = gca;
% ax.Position
% ax.XLim
% 
% position_aux = ax.Position
% ratio = ax.XLim(2)/ax.YLim(2)
% ax.Position(3) = ax.Position(4) * ratio
% 
% RX_x = [RX_x_2E_open_1;RX_x_2E_open_2;RX_x_corridor;RX_x_2S_open_1;RX_x_2S_open_2];
% RX_y = [RX_y_2E_open_1;RX_y_2E_open_2;RX_y_corridor;RX_y_2S_open_1;RX_y_2S_open_2];
% 
% points_x = [points_x_2E_open;points_x_marconi;points_x_2E_corridor.';points_x_corridor.';point_x_area.';points_x_2S_open;points_x_2E_office];
% points_y = [points_y_2E_open;points_y_marconi;points_y_2E_corridor.';points_y_corridor.';point_y_area.';points_y_2S_open;points_y_2E_office];
% 
% RX_x = RX_x([4,3,2,5,1]);
% RX_y = RX_y([4,3,2,5,1]);
% 
% fig_2 = figure('units','normalized','outerposition',[0 0 0.4 0.8])
% prx1 = plot(RX_x,RX_y, "r*", "LineWidth", 4)
% hold on
% p_points = plot(points_x,points_y, "b*", "LineWidth", 4)
% 
% text(points_x + 0.2, points_y + 0.2,string(1:length(points_x)))
% text(RX_x + 0.2, RX_y + 0.2,string([4,5,6,7,10]))
% 
% distances_RX = sqrt((points_x - RX_x.').^2 + (points_y - RX_y.').^2);
% diff_RX_x = points_x - RX_x.';
% diff_RX_y = points_y - RX_y.';
% angles_RX = rad2deg(atan(diff_RX_x./diff_RX_y));
% angles_RX_2 = rad2deg(atan(diff_RX_y./diff_RX_x));
% 
% angles_RX(:,[3 5]) = angles_RX_2(:,[3 5]);
% angles_RX(:,2) = rad2deg(atan((-1)*diff_RX_y(:,2)./diff_RX_x(:,2)));
% angles_RX(45:2:57,2) = angles_RX(45:2:57,2)*(-1);
% 
% mkdir("../../mat_files/Data_True");
% 
% save("../../mat_files/Data_True/data_distance_angle_true", "angles_RX", "distances_RX")
% save("../../mat_files/Data_True/points_coordinates", "points_x", "points_y")
% save("../../mat_files/Data_True/RX_coordinates", "RX_x", "RX_y")
% 
% % fig_2.CurrentAxes
% % ax = gca;
% % ax.Position
% % ax.XLim
% % 
% % position_aux = ax.Position
% % ratio = ax.XLim(2)/ax.YLim(2)
% % ax.Position(3) = ax.Position(4) * ratio
% 
% print(fig, "Map_Floop_Plan", "-dpng")
% 
% area_m2 = distance_y_2S_corridor * (limit_x_2S_open + offices_2S) + (limit_y_corridor * limit_x_corridor) + (limit_y_2E_open * (limit_x_2E_corridor + limit_x_2E_marconi*2))