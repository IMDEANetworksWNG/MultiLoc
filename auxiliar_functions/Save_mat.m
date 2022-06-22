%% M-Cube antennas

antenna_positions_mcube = [nan, 24, 22, 6, 8, nan; 19, 17, 23, 7, 1, 3; 20, 18, 21, 5, 2, 4; 28, 27, 26, 9, 11, 12; 25, 29, 32, 16, 13, 10; nan, 31, 30, 14, 15, nan];
antenna_positions       = rot90(33-antenna_positions_mcube);
%% Process the data
out_path = 'processed_data/Measurements_IMDEA/';
out_name = 'pan_0_tilt_0';
file     = 'raw_data/grid_data/pan_0_tilt_0.txt';

[magnitudes, phases, times] = Parse_raw_data(file);

%% Save to a .mat
save([ out_path out_name '.mat'], 'magnitudes', 'phases', 'times', 'antenna_positions')
