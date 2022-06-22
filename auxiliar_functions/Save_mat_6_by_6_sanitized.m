clc 
clear all
close all

%% M-Cube antennas

% Load antenna data
load('/home/imdea/Documents/MATLAB/Music_Mikrotik/processed_data/antennas_mikrotik.mat')
%% Process the data
out_path = 'processed_data/oscillator_calibration_measurements/44/';
in_path  = 'raw_data/oscillator_calibration_measurements/44/';

for pan=-10:5:10
    for tilt=-10:5:10
        
        name = ['pan_' num2str(pan) '_tilt_' num2str(tilt)];
        
        file = [in_path name  '.txt'];

        [magnitudes, phases, times] = Parse_raw_data(file);

        % Transform to a 6x6 matrix
        num_samples = 300;

        % Clean the data
        pre_channel    = zeros(6, 6, num_samples);
        pre_magnitudes = zeros(6, 6, num_samples);

        % We go up to 30 instead of 32
        % so that 31 and 32 are disabled
        % since they return random data
        for jj=1:30

            a = phases(:, jj);
            a = a*2*pi/1024;
            % move to complex
            a = exp(1i*a);

            try
                [a, phase_offset_0] = Sanitize(a);
            catch
                disp(['Converging error on file ' file])
                broken = 1;
                break
            end

            % Remove oscilator
            a = a/exp(1i*antenna_oscilator_phases(antenna_positions == jj));

            [row,col] = find(antenna_positions == jj);
            pre_channel(row, col, :) = a(1:num_samples, :);

            % Magnitudes
            b = magnitudes(:, jj);

            pre_magnitudes(row, col, :) = b(1:num_samples, :);
        end

        magnitudes = pre_magnitudes;
        phases     = pre_channel;
        % Save to a .mat
        save([ out_path name '.mat'], 'magnitudes', 'phases', 'times', 'antenna_positions')

    end
end
