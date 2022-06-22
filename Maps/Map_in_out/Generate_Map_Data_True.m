clear
clc
close all

RX_x = [0];
RX_y = [0];

dis_in_out = 0.14 + 1.45 + 1.34;
dis_in_out2 = 0.3*9 + 0.64;

points_x = reshape(repmat(-2:2,4,1).',[],1)*2.97;
points_y = [ones(5,1)*4,ones(5,1)*(8), ones(5,1) *(8 + dis_in_out), ones(5,1) *(8 + dis_in_out+dis_in_out2)];
points_y = points_y(:);
points_x(15) = points_x(15) + 0.7;

RX_x(2) = points_x(end) + 0.16 + 0.13 + 16*0.4;
RX_y(2) = points_y(end) - (0.08 + 4*0.4);

diff_x = RX_x - points_x;
diff_y = RX_y - points_y;

distances_RX = sqrt(diff_x.^2 + diff_y.^2);
angle = angle(diff_x + 1i*diff_y);
angles_RX = rad2deg(angle);
angles_RX(:,1) = angles_RX(:,1) + 90;

figure, plot(RX_x, RX_y, "*r")
hold on
plot(points_x, points_y, "*b")

name_points = string(1:20)
text(points_x + 0.2, points_y + 0.2,name_points)

name_RX = ["10/11/44"; "4/12/43"];
text(RX_x + 0.2, RX_y + 0.2,name_RX)


mkdir("../../mat_files")
mkdir("../../mat_files/in_out_map")


save("../../mat_files/in_out_map/data_distance_angle_true", "distances_RX", "angles_RX")
save("../../mat_files/in_out_map/points_coordinates", "points_x", "points_y")
save("../../mat_files/in_out_map/RX_coordinates", "RX_x", "RX_y")

savefig("../../mat_files/in_out_map/map")
