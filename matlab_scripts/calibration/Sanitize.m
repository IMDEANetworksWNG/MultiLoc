function [csi_data,offset,converged] = Sanitize(csi_data)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% take the data for antenna 2
%     csi_data = csi_data(:,17);

    % create a grid for the possible values
    grid = linspace(-pi,pi,1024);

    % move to complex number
    S_phase = exp(1i*grid);

    % create the spectrum. To do that, we multiply by the conjugate. The phase
    % that maximizes the real part is the correct value
    spectrum_phase = csi_data * conj(S_phase);

    % take the max
    [~, index_max] = max(sum(real(spectrum_phase),1));

    % we center based on the angle that maximizes the real part
    offset_phase = grid(index_max);

    csi_data = csi_data .* conj(S_phase(index_max));

    csi_data = csi_data .* exp(1i*pi/2);


    % apply guaxian mixture model
    GM = fitgmdist(angle(csi_data),2);

    % show the mean of the gaussian distribution
    GM.mu;
    GM.ComponentProportion;
    converged = GM.Converged;
    % order the phase in term of the proportion
    [~, index_sort] = sort(GM.ComponentProportion, "descend");

    mu_sorted = GM.mu(index_sort);
    % take the phase in the middle
    phase_middle = mean(mu_sorted);

    % sum to the outliers pi/2
    if (mu_sorted(2) < mu_sorted(1))
        index_out = angle(csi_data) < phase_middle;
    else
        index_out = angle(csi_data) > phase_middle;
    end

    % check whether is pi/2, pi and 3pi/2
    grid_pi = [pi/2, pi, -pi/2];

    diff_mu = diff(mu_sorted);

    % look for the minimum of the substraction
    sub_diff_mu = diff_mu + grid_pi;
    [~, index_min_mu] = min(abs(sub_diff_mu));


    % analizy when the jump is pi/2, pi or 3*pi/2
    csi_data(index_out) = csi_data(index_out) .* exp(1i*grid_pi(index_min_mu));

    % csi_data = csi_data .* exp(-1i*offset_mu);


    % remove the pi/2 offset
    csi_data = csi_data .* exp(-1i*pi/2);

    % it should be center at 0 or closer to 0
    mean_offset = median(angle(csi_data));

    csi_data = csi_data .* exp(-1i*mean_offset);

    offset = mean_offset + offset_phase;

    csi_data = csi_data*exp(1i*offset);
end
