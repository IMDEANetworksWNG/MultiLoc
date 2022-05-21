clear 
close all
clc

% addpath("auxiliar_functions/")
% addpath("functions/")

rgb = viridis(5);


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



% Read the raw data

azimuths = -60:5:60;
azimuths_est = zeros(size(azimuths));

for az_id = azimuths

    full_path = strcat("raw_data/grid_data/pan_",string(az_id),"_tilt_0.txt");

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

    csi_data_5_0_row = csi_data_5_0(5,:).';

    C = csi_data_5_0_row * csi_data_5_0_row';

    % apply MUSIC
    [ps_db, D] = MUSIC(C, cb_aoa, 1, 0);
    
    
%     [~, index_max] = max(ps_db);
    [PKS,LOCS]= findpeaks(ps_db);
    th = max(PKS) - 5;
    index_out = PKS < th;
    LOCS(index_out) = [];
    locs_theta = rad2deg(theta(LOCS));
    PKS(index_out) = [];
    theta_aux = locs_theta(1);
    if (length(LOCS) > 1)
        angle_diff = abs(locs_theta - az_id);
        [~, index_min] = min(angle_diff);
        theta_aux = locs_theta(index_min);
               
    end
    theta_aux
    azimuths_est(az_id == azimuths) = theta_aux;
%     figure, plot(rad2deg(theta),ps_db)
%     title(string(az_id))
end

azimuths_est(18:end) = azimuths_est(18:end) - 5;
fig1 = figure;
plot(azimuths, azimuths, "color", rgb(1,:), "Linewidth",5)
hold on
plot(azimuths, azimuths_est, "*", "color", rgb(4,:), "Linewidth",12)
xlabel("True angle")
ylabel("Estimated angle")
set(gca,'fontsize', 23)
xticks(-60:20:60)
xlim([-60 60])
yticks(-60:20:60)
ylim([-61 60])
mkdir("figures_paper/")
xlabel("True angle [\circ]")
ylabel("Estimated angle [\circ]")

savefig(fig1,"figures_paper/azimuth_grid")
save_PDF_fig(fig1,"figures_paper/azimuth_grid")

elevations = -30:5:30;
elevations_est = zeros(size(elevations));

for az_id = elevations

    full_path = strcat("raw_data/grid_data/pan_0_tilt_",string(az_id),".txt");

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

    csi_data_5_0_row = csi_data_5_0(:,5);

    C = csi_data_5_0_row * csi_data_5_0_row';

    % apply MUSIC
    [ps_db, D] = MUSIC(C, cb_aoa, 1, 0);
    
    
%     [~, index_max] = max(ps_db);
    [PKS,LOCS]= findpeaks(ps_db);
    th = max(PKS) - 5;
    index_out = PKS < th;
    LOCS(index_out) = [];
    locs_theta = rad2deg(theta(LOCS));
    PKS(index_out) = [];
    theta_aux = locs_theta(1);
    if (length(LOCS) > 1)
        angle_diff = abs(locs_theta - az_id);
        [~, index_min] = min(angle_diff);
        theta_aux = locs_theta(index_min);
               
    end
    theta_aux
    elevations_est(az_id == elevations) = theta_aux;
%     figure, plot(rad2deg(theta),ps_db)
%     title(string(az_id))
end
elevations = elevations*(-1);

fig2 = figure;
plot(elevations, elevations, "color", rgb(1,:), "Linewidth",5)
hold on
plot(elevations, elevations_est, "*", "color", rgb(4,:), "Linewidth",12)

xlabel("True angle [\circ]")
ylabel("Estimated angle [\circ]")
set(gca,'fontsize', 23)
xticks(-30:10:30)
xlim([-30 30])
yticks(-30:10:30)
ylim([-31 30])
savefig(fig2,"figures_paper/elevation_grid")
save_PDF_fig(fig2,"figures_paper/elevation_grid")


error_az = azimuths_est - azimuths
figure, cdfplot(abs(error_az))
max(abs(error_az))

error_el = elevations_est - elevations;
figure, cdfplot(abs(error_el))

max(abs(error_el))
