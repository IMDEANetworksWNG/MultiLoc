function [index_az_max,index_el_max,attenuation] = Individual_Azimuth_Elevation_Estimator(channel,S_az,S_el,index_el)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
        
    % estimate azimuth knowing elevation
    channel_el = sum(channel .* S_el(:,index_el)',2);
    spectrum_az = (S_az' * channel_el);
    [~,index_az_max] = max(abs(spectrum_az));

    % estimate elevation knowing azimuth
    channel_az = sum(channel .* conj(S_az(:,index_az_max)),1);
    spectrum_el = (S_el' * channel_az.');
    [~,index_el_max] = max(abs(spectrum_el));
    attenuation = spectrum_el(index_el_max);

    
end

