function [Az_estimated, El_estimated, att] = mD_track_2D(channel, S_az, S_el)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    % number of antennas and active subcarriers
    [N, M] = size(channel);
    
    
    index_zeros = channel == 0;

    %% In mD-track paper is 3.2.1 Initial estimation

    % take auxiliar variables to remove them in the loop
    channel_aux = channel;

    % variables for the arguments of angles and times of arrivas
    arg_max_el = zeros(1,0);
    arg_max_az = zeros(1,0);

    % variable for the power of the attenuation
    power = zeros(1,0);
    att = zeros(1,0);

      
    % variable for reconstract the estimate parameters
    channel_recontract = zeros(0,N,M);

    out = false;
    l = 0;
    while (out == 0)
        % iterate over the number of paths
        l = l + 1;
        % z functions
        [matrix_az_el] = Jointly_Azimuth_Elevation_Estimator(S_az, S_el, channel_aux);
        
        % take the power
        matrix_az_el_power = abs(matrix_az_el).^2;
        max_value = max(matrix_az_el_power(:));

        % take the position of the maximum
        [index_az_max, index_el_max] = find(matrix_az_el_power == max_value);

        if length(index_az_max) > 1
            index_az_max = index_az_max(1);
            index_el_max = index_el_max(1);
        end
        
        attenuation = matrix_az_el(index_az_max, index_el_max);
        attenuation = attenuation/((N*M)-sum(index_zeros(:)));
        att(l) = attenuation;

        % take the power of the maximum
        power(l) = abs(attenuation).^2;
        
        
        % check the power
        if (power(1) * 10^(-25/10) > power(l)) % Originally (-25/10)
            out = true;
            power(l) = [];
            att(l) = [];
            l = l -1;
            break;
        end
        
        % take the argument of the path to be saved
        arg_max_el(l) = index_el_max;
        arg_max_az(l) = index_az_max;

        % take the pattern of the component
        pattern_el = (S_el(:,arg_max_el(l)));
        pattern_az = S_az(:,arg_max_az(l));
       
        % create the component
        channel_remove_el = repmat(pattern_el,1,N).';
%         channel_remove_toa = channel_remove_toa * attenuation;
        channel_remove_el_az = (channel_remove_el .* pattern_az)*attenuation;

        channel_remove_el_az(index_zeros) = 0;
        % remove it from the channel
        channel_aux = channel_aux - channel_remove_el_az;

        % recontract the signal
        channel_recontract(l,:,:) = channel_remove_el_az;
    end
    
   
    % number of estimated paths
    L_estimated = l;
    
    %% In mD-track paper is 3.2.2  Iterative path parameter refinement.

    % take the output of the previous part as a residual
    channel_residual = channel_aux;

    % fix to ten iterations
    for iteration = 1:10
    %     iteration
        for l = 1:(L_estimated)

            % add noise
            channel_component = squeeze(channel_recontract(l,:,:)) + channel_residual;

            [index_az_max,index_el_max,attenuation] = Individual_Azimuth_Elevation_Estimator(channel_component,S_az,S_el,arg_max_el(l));
            attenuation = attenuation/((N*M)-sum(index_zeros(:)));
            att(l) = attenuation;

            % take the power
            power(l) = abs(attenuation).^2;
            
            arg_max_el(l) = index_el_max;
            arg_max_az(l) = index_az_max;

            % take the pattern of the component
            pattern_el = (S_el(:,arg_max_el(l)));
            pattern_az = S_az(:,arg_max_az(l));
                       
            % reconstract the signal
            channel_remove_el = repmat(pattern_el,1,N).';
%             channel_remove_toa = channel_remove_toa .* (cos(shift_phase_component(l)) + 1i*sin(shift_phase_component(l)));
%             channel_remove_toa = channel_remove_toa * attenuation;

            channel_remove_el_az = channel_remove_el .* pattern_az;
            channel_remove_el_az(index_zeros) = 0;
            channel_recontract(l,:,:) = channel_remove_el_az * attenuation;

%             channel_recontract(l,:,:) = channel_remove_toa_aoa * sqrt(power(l));

            
            % remove everything from the channel
            channel_sum = sum(channel_recontract,1);
            channel_sum = squeeze(channel_sum);
            channel_residual = channel - channel_sum;

        end
        
    end
    
    
    % take toa and aoa
%     ToA_estimated = times_toa(arg_max_toas);
%     AoA_estimated = rad2deg(angles(arg_max_aoas));
    El_estimated = arg_max_el;
    Az_estimated = arg_max_az;
end

