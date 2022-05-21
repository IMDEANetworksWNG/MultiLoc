function [matrix_az_el] = Jointly_Azimuth_Elevation_Estimator(S_az, S_el, channel_aux)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [~,len_az] = size(S_az);
    [~,len_el] = size(S_el);
    matrix_az_el = zeros(len_az, len_el);
    
    
    for index_az = 1:len_el
        
        channel_az = S_az(:,index_az)'.*channel_aux.';
        channel_aoa = sum(channel_az,2);

%             channel_freq = fft(channel_aoa);
        channel_az_el = (S_el')*channel_aoa;
        alpha = sum((channel_az_el),2);

        matrix_az_el(index_az, :) = alpha;
    end
end

