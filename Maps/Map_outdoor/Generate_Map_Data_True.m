close all
clc
clear

%% Points
tile_size = 0.4; % In meters

labels = [111:126];

points_x = [10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 0, 15, 30, -6, -6, -6];
points_y = [0, -10, -20, -30, -40, 0, -10, -20, -30, -40, 8, 8, 8, 0, -15, -30];

% AP_labels = [43:46];
AP_labels_LF = [4, 6, 7, 10];
AP_labels_HF = [43, 44, 45, 46];

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

text(RX_x + 0.1, RX_y + 0.1, strcat(string(AP_labels_LF), "/", string(AP_labels_HF)));

p_points = plot(points_x,points_y, "b*", "LineWidth", 4);
text(points_x + 0.1, points_y + 0.1, string(labels));

% take the distances between APs to points
distances_RX = sqrt((points_x' - RX_x).^2 + (points_y' - RX_y).^2);

% take the angles between APs to points
diff_RX_x = points_x' - RX_x;
diff_RX_y = points_y' - RX_y;

points_complex = diff_RX_x + 1i*diff_RX_y;

angles_RX = angle(points_complex);
angles_RX(angles_RX < 0) = angles_RX(angles_RX < 0) + 2*pi;
angles_RX = rad2deg(angles_RX);

% move 90 to 0
angles_RX = angles_RX - 90;

% apply the rotation
rotation = [-45 -135 -225 -90]

angles_without_rotations = angles_RX;

angles_RX = angles_RX + rotation;
angles_RX = angles_RX(:, [1, 2, 3, 4]);
distances_RX = distances_RX(:, [1, 2, 3, 4]);

% Point labels
point_labels = 111:126;

mkdir("../../mat_files")
mkdir("../../mat_files/outdoor_map")

save("../../mat_files/outdoor_map/points_coordinates", "points_x", "points_y", "point_labels")
save("../../mat_files/outdoor_map/RX_coordinates", "RX_x","RX_y")
save("../../mat_files/outdoor_map/data_distance_angle_true", "angles_RX", "distances_RX", "AP_labels_LF", "AP_labels_HF")

%% 
% Calculate the optimal rotation for HF
rotations = sort([0, 45, 270, 315, 90, 135, 180, 225]);

% We map it from [-90, 90] to [0, 360], we multiply by -1
% so it is from the POV of the point
outdoor_true_angles_360 = mod(-1*angles_without_rotations, 360);

% Since the routers are pointing down we need to rotate the angles
outdoor_true_angles_360 = mod(360-outdoor_true_angles_360, 360);

% Now we need to round to the closest multiple of 45
outdoor_true_angles_360 = round(outdoor_true_angles_360 / 45) * 45;

% And 360 is equal to 0
outdoor_true_angles_360(outdoor_true_angles_360==360) = 0;

% Convert the angles to indexes in the matrix
for ii=1:numel(outdoor_true_angles_360)
   
    outdoor_true_angles_360(ii) = find(rotations==outdoor_true_angles_360(ii));
end

optimal_rotation_index = outdoor_true_angles_360;