function [ToA_estimated, AoA_estimated, AoD_estimated, power] = mD_track_3D(channel, step_toa, step_aoa, step_aod, K, N, M, BW, grid_toa)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % take the steering vectors in time and in angle domain
    % S_toa represents the phase difference due to a delay
    % S_aoa represents the phase difference due to an angle of arrival
    % S_aod represents the phase difference due to an angle of departure
    % times and angles_aoa angles_aod are the grid for times and angles of
    % arrival and departure
    
    [S_toa, S_aoa, S_aod, times, angles_aoa, angles_aod] = Grid_ToA_AoA_AoD_Grid(step_toa,step_aoa, step_aod, K, N, M, BW, grid_toa);

    % do a rough estimation of time of arrival to chech whether is the
    % power in time
    channel_toa = squeeze(channel(:,:,1));
   
    % compute the power
    channel_steering = (S_toa')*channel_toa;   
    ps_db_toa = sum(abs(channel_steering),2);

    % take the maximum
    [~, LOCKS] = max(ps_db_toa);
    
    % take between the maximum a window of 16 samples
    index_toa = zeros(16/step_toa,1);
    index_toa(:,1) = ((LOCKS-1)- (8/(step_toa))):((LOCKS-2) + (8/(step_toa)));
    index_toa(:,1) = mod(index_toa, length(S_toa))+1;
    
    % recompute S_toa
    S_toa = S_toa(:,index_toa(:));
    times_toa = times(index_toa(:));
    %% In mD-track paper is 3.2.1 Initial estimation

    % take auxiliar variables to remove them in the loop
    channel_aux = channel;

    % variables for the arguments of angles of arrival and departure and times of arrivas
    arg_max_toas = zeros(1,0);
    arg_max_aoas = zeros(1,0);
    arg_max_aods = zeros(1,0);

    % variable for the power of the attenuation
    power = zeros(1,0);

    % phase of the attenuation
    shift_phase_component = zeros(1,0);
    
    % number of active subcarriers
    K_active = length(grid_toa);
    
    % variables for reconstract the signal
    channel_recontract = zeros(0,K_active, N, M);

    out = false;
    l = 0;
    while (out == 0)   
        l = l + 1;
        % estimate the power for all the toas and aoas
        [matrix_toa_aoa_aod] = Jointly_ToA_AoA_AoD_Estimator_3(S_toa,S_aoa,S_aod, channel_aux);

%         [matrix_toa_aoa_aod_2] = Jointly_ToA_AoA_AoD_Estimator_3(S_toa,S_aoa,S_aod, channel_aux);
        
        matrix_toa_aoa_aod_power = abs(matrix_toa_aoa_aod).^2;
        % take into account if the peaks in the power are at the end and at
        % the beginning of the time
        change_Data = 0;
        lim_down = times_toa(1);
        lim_up = times_toa(end);
        
        if(lim_down > lim_up)
            
            times_aux = circshift(times, -length(index_toa));
            times_toa = times_aux(index_toa(:));
            lim_down = times_toa(1);
            lim_up = times_toa(end);
        end
        
        % take the power
        [power(l), index_max] = max(matrix_toa_aoa_aod_power(:));
        
        % check the power
        if (power(1) * 10^(-25/10) > power(l) || l > 25)
            out = true;
            power(l) = [];
            l = l - 1;
            break;
        end
        
        [index_toa_max, index_aoa_max, index_aod_max] = ind2sub(size(matrix_toa_aoa_aod_power),index_max);

        if length(index_aoa_max) > 1
            index_aoa_max = index_aoa_max(1);
            index_aod_max = index_aod_max(1);
            index_toa_max = index_toa_max(1);
        end
        
        arg_max_aods(l) = index_aod_max;
        arg_max_toas(l) = index_toa_max;
        arg_max_aoas(l) = index_aoa_max;

        % take the pattern of the component
        pattern_toa = (S_toa(:,arg_max_toas(l)));
        pattern_aoa = S_aoa(:,arg_max_aoas(l));
        pattern_aod = (S_aod(:,arg_max_aods(l)));

        % estimate the phase of the attenuation by iteration 
        shift_phase_component(l) = angle(matrix_toa_aoa_aod(index_toa_max, index_aoa_max, index_aod_max));
        
        % create the component
        channel_remove_toa = repmat(pattern_toa,1,N,M);
        channel_remove_toa = channel_remove_toa .* (cos(shift_phase_component(l)) + 1i*sin(shift_phase_component(l)));
        channel_remove_toa_aoa = channel_remove_toa .* pattern_aoa.';
        channel_remove_toa_aoa_aod = channel_remove_toa_aoa .* reshape(pattern_aod.',1,1,M);



        % remove it from the channel
        channel_remove_toa_aoa_aod_power = channel_remove_toa_aoa_aod * sqrt(power(l));
        channel_aux = channel_aux - channel_remove_toa_aoa_aod_power;

        % recontract the first signal from the signal
        channel_recontract(l,:,:,:) = channel_remove_toa_aoa_aod_power;
        


    end
    
    L_estimated = l;
    %% In mD-track paper is 3.2.2  Iterative path parameter refinement.

    % take the output of the previous part as a residual


    channel_residual = channel_aux;


    % fix to ten iterations

    for repetition = 1:10
    %     l
        for l = 1:(L_estimated)

            % add noise
            channel_component = squeeze(channel_recontract(l,:,:,:)) + channel_residual;
            
            % estimate the power for all the toas and aoas and aods
            % estimate aoa fixing toa and aod
            [alpha_aoa] = Jointly_ToA_AoA_AoD_Estimator_2(S_toa(:,arg_max_toas(l)),S_aoa,S_aod(:,arg_max_aods(l)), channel_component);
            [~, index_aoa_max] = max(abs(alpha_aoa));

            % estimate toa fixing aoa and aod
            [alpha_toa] = Jointly_ToA_AoA_AoD_Estimator_2(S_toa,S_aoa(:,index_aoa_max),S_aod(:,arg_max_aods(l)), channel_component);
            [~, index_toa_max] = max(abs(alpha_toa));

            % estimate aod fixing aoa and toa
            [alpha_aod] = Jointly_ToA_AoA_AoD_Estimator_2(S_toa(:,index_toa_max),S_aoa(:,index_aoa_max),S_aod, channel_component);
            [~, index_aod_max] = max(abs(alpha_aod));

            
            [attenuation] = Jointly_ToA_AoA_AoD_Estimator_2(S_toa(:,index_toa_max),S_aoa(:,index_aoa_max),S_aod(:,index_aod_max), channel_component);

            power(l) = abs(attenuation).^2;
            

            if length(index_aoa_max) > 1
                index_aoa_max = index_aoa_max(1);
                index_aod_max = index_aod_max(1);
                index_toa_max = index_toa_max(1);
            end

            arg_max_aods(l) = index_aod_max;
            arg_max_toas(l) = index_toa_max;
            arg_max_aoas(l) = index_aoa_max;

            % take the pattern of the component
            pattern_toa = (S_toa(:,arg_max_toas(l)));
            pattern_aoa = S_aoa(:,arg_max_aoas(l));
            pattern_aod = (S_aod(:,arg_max_aods(l)));
            
            shift_phase_component(l) = angle(attenuation);
            

            
            % create the component
            channel_remove_toa = repmat(pattern_toa,1,N,M);
            channel_remove_toa = channel_remove_toa .* (cos(shift_phase_component(l)) + 1i*sin(shift_phase_component(l)));
            channel_remove_toa_aoa = channel_remove_toa .* pattern_aoa.';
            channel_remove_toa_aoa_aod = channel_remove_toa_aoa .* reshape(pattern_aod.',1,1,M);


            channel_remove_toa_aoa_aod_power = channel_remove_toa_aoa_aod * sqrt(power(l));
            channel_recontract(l,:,:,:) = channel_remove_toa_aoa_aod_power;
            
            % remove everything from the channel
            channel_sum = sum(channel_recontract,1);
            channel_sum = squeeze(channel_sum);
            channel_residual = channel - channel_sum;
        end
        
    end
    
    

    ToA_estimated = times_toa(arg_max_toas);
    AoA_estimated = rad2deg(angles_aoa(arg_max_aoas));
    AoD_estimated = rad2deg(angles_aod(arg_max_aods));   
end

