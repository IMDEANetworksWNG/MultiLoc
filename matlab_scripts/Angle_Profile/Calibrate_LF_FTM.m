load('mat_files/FTM/LF/estimated_distance.mat');
load('mat_files/indoor_map/data_distance_angle_true.mat');

error = sqrt((distances_RX - estimated_distance).^2);

cdfplot(error(:))
