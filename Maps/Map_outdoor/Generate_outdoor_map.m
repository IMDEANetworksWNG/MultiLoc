%% Points
tile_size = 0.4; % In meters

labels = [111:126];

points_x = [10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 0, 15, 30, -6, -6, -6];
points_y = [0, -10, -20, -30, -40, 0, -10, -20, -30, -40, 8, 8, 8, 0, -15, -30];

AP_labels = [43:46];

RX_x = [39, 34, -5, 39];
RX_y = [-46, 15, 15, -23];

% We convert everything to meters
points_x = points_x * tile_size;
points_y = points_y * tile_size;

RX_x = RX_x * tile_size;
RX_y = RX_y * tile_size;

% Aps 43 and 46 are 9 cm off on the X axis, fix that
RX_x(1) = RX_x(1) - 0.09;
RX_x(4) = RX_x(4) - 0.09;

% Move everything so that (0, 0) is bottom left
y_offset = abs(RX_y(1));
x_offset = abs(points_x(16));
RX_y = RX_y + y_offset;
RX_x = RX_x + x_offset;

points_x = points_x + x_offset;
points_y = points_y + y_offset;

%% Generate the map

figure,

% Grass
rectangle('Position', [0 + x_offset -50*tile_size+y_offset 34*tile_size 50*tile_size], 'EdgeColor','g');
hold on

p_points = plot(points_x,points_y, "b*", "LineWidth", 4);
p_routers = plot(RX_x,RX_y, "r*", "LineWidth", 4);

text(RX_x + 0.1, RX_y + 0.1, string(AP_labels));

p_points = plot(points_x,points_y, "b*", "LineWidth", 4);
text(points_x + 0.1, points_y + 0.1, string(labels));

% take the distances between APs to points
true_distances = sqrt((points_x' - RX_x).^2 + (points_y' - RX_y).^2);

% take the angles between APs to points
diff_RX_x = points_x' - RX_x;
diff_RX_y = points_y' - RX_y;
true_angles = rad2deg(atan(diff_RX_x./diff_RX_y));
