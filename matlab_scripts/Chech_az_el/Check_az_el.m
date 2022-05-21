clear 
close all
clc

% frequency
freq = 60.48e9;

% speed of light
c = 3e8;

% the wavelength
lambda = c/freq;

% distance between antennas
d = lambda*0.58;

step_angle = 0.5;

N = 6;
% create the codebo3ok
[cb_aoa, theta] = Grid_AoA(step_angle, N,d,lambda);

% Load antenna data
load('processed_data/antennas_mikrotik.mat')


% Read the raw data
full_path = ['raw_data/grid_data/pan_0_tilt_0.txt'];

[magnitudes, phases, ~] = Parse_raw_data(full_path);

broken = 0;

[num_samples,~] = size(magnitudes);

% Clean the data
pre_channel = zeros(6, 6, num_samples);

% We go up to 30 instead of 32
% so that 31 and 32 are disabled
% since they return random data
for jj=1:30

    a = phases(:, jj);
    a = a*2*pi/1024;
    % move to complex
    a = exp(1i*a);

    converging_limit = 50;
    converging_retries = 0;
    converged = 0;
    while converged == 0

        try
            [a, phase_offset_0, converged] = Sanitize(a);
%                         figure; plot(angle(a));
            %broken = 0;
        catch
            disp(['Converging error on file ' full_path])
            %broken = 1;
            %break
        end

        if converging_retries == converging_limit

            break
        end

        converging_retries = converging_retries + 1;

    end

    if converging_retries == converging_limit
        disp(['Converging threshold reached, ignoring ' filename])
        continue
    end
    
    [row,col] = find(antenna_positions == jj);
    pre_channel(row, col, :) = a(1:num_samples, :);
end

csi_data = pre_channel;
csi_data_osc = sum(csi_data,3)/num_samples;

csi_data_0_0 = csi_data./csi_data_osc;


% Read the raw data
full_path = ['raw_data/grid_data/pan_30_tilt_20.txt'];

[magnitudes, phases, ~] = Parse_raw_data(full_path);

broken = 0;

[num_samples,~] = size(magnitudes);

% Clean the data
pre_channel = zeros(6, 6, num_samples);

% We go up to 30 instead of 32
% so that 31 and 32 are disabled
% since they return random data
for jj=1:30

    a = phases(:, jj);
    a = a*2*pi/1024;
    % move to complex
    a = exp(1i*a);

    converging_limit = 50;
    converging_retries = 0;
    converged = 0;
    while converged == 0

        try
            [a, phase_offset_0, converged] = Sanitize(a);
%                         figure; plot(angle(a));
            %broken = 0;
        catch
            disp(['Converging error on file ' full_path])
            %broken = 1;
            %break
        end

        if converging_retries == converging_limit

            break
        end

        converging_retries = converging_retries + 1;

    end

    if converging_retries == converging_limit
        disp(['Converging threshold reached, ignoring ' filename])
        continue
    end
    
    [row,col] = find(antenna_positions == jj);
    pre_channel(row, col, :) = a(1:num_samples, :);
end

csi_data = pre_channel;
csi_data = sum(csi_data,3)/num_samples;

csi_data_5_0 = csi_data./csi_data_osc;

csi_data_5_0_row = csi_data_5_0(3,:).';

C = csi_data_5_0_row * csi_data_5_0_row';

% apply MUSIC
[ps_db, D] = MUSIC(C, cb_aoa, 1, 0);

figure, plot(rad2deg(theta),ps_db)

for ii = 1:N
    cb_aux(ii,:) = exp(-1i*(ii-1)*d*2*pi*sin(theta)*cosd(-20)/lambda);
end

% apply MUSIC
[ps_db, D] = MUSIC(C, cb_aux, 1, 0);

figure, plot(rad2deg(theta),ps_db)

angles = unwrap(angle(csi_data_5_0_row./csi_data_5_0_row(1)))


% Read the raw data
full_path = ['raw_data/grid_data/pan_-30_tilt_25.txt'];

[magnitudes, phases, ~] = Parse_raw_data(full_path);

broken = 0;

[num_samples,~] = size(magnitudes);

% Clean the data
pre_channel = zeros(6, 6, num_samples);

% We go up to 30 instead of 32
% so that 31 and 32 are disabled
% since they return random data
for jj=1:30

    a = phases(:, jj);
    a = a*2*pi/1024;
    % move to complex
    a = exp(1i*a);

    converging_limit = 50;
    converging_retries = 0;
    converged = 0;
    while converged == 0

        try
            [a, phase_offset_0, converged] = Sanitize(a);
%                         figure; plot(angle(a));
            %broken = 0;
        catch
            disp(['Converging error on file ' full_path])
            %broken = 1;
            %break
        end

        if converging_retries == converging_limit

            break
        end

        converging_retries = converging_retries + 1;

    end

    if converging_retries == converging_limit
        disp(['Converging threshold reached, ignoring ' filename])
        continue
    end
    
    [row,col] = find(antenna_positions == jj);
    pre_channel(row, col, :) = a(1:num_samples, :);
end

csi_data = pre_channel;
csi_data = sum(csi_data,3)/num_samples;

csi_data_5_0 = csi_data./csi_data_osc;

csi_data_5_0_row = csi_data_5_0(:,4);

C = csi_data_5_0_row * csi_data_5_0_row';

% apply MUSIC
[ps_db, D] = MUSIC(C, cb_aoa, 1, 0);

figure, plot(rad2deg(theta),ps_db)

angles = unwrap(angle(csi_data_5_0_row./csi_data_5_0_row(1)));

for ii = 1:N
    cb_aux(ii,:) = exp(-1i*(ii-1)*d*2*pi*sin(theta)*cosd(25)/lambda);
end

% apply MUSIC
[ps_db, D] = MUSIC(C, cb_aux, 1, 0);

figure, plot(rad2deg(theta),ps_db)

% unwrap(angle(cb_aux)).'