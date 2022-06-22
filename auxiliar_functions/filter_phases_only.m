%file = 'raw_data/with_rotation_table/pan_4.0_tilt_0.txt';

function [offset_phase, csi_data_test_aux] = filter_phases_only(csi_data_test)

    % move the phase to radians
    %phases = (phases/1024)*2*pi;

    % move to complex
    %csi_data = exp(1i*phases);


    % take the data for antenna 2
    %csi_data_test = csi_data(:, antenna);
    %figure, plot(angle(csi_data_test));
    %figure, cdfplot(angle(csi_data_test));
    %figure, histogram(angle(csi_data_test), 100, "Normalization", "probability");

    % take the mode and the interqualtic mean
    mode_phase = mode(angle(csi_data_test));


    % create a grid for the possible values
    grid = linspace(-pi,pi,1024);

    % move to complex number
    S_phase = exp(1i*grid);
    % plot it, the phase should be linar
    %figure, plot(angle(S_phase))

    % create the spectrum. To do that, we multiply by the conjugate. The phase
    % that maximizes the real part is the correct value
    spectrum_phase = csi_data_test * conj(S_phase);
    %figure, plot(grid,sum(real(spectrum_phase),1))

    % take the max
    [~, index_max] = max(sum(real(spectrum_phase),1));
    grid(index_max);

    % we center based on the angle that maximizes the real part
    offset_phase = grid(index_max);

    csi_data_test_aux = csi_data_test .* conj(S_phase(index_max));
    %figure, plot(angle(csi_data_test_aux));

    csi_data_test_aux = csi_data_test_aux .* exp(1i*pi/2);


    %figure, histogram(angle(csi_data_test_aux), 100, "Normalization", "probability");


    % apply guaxian mixture model
    GM = fitgmdist(angle(csi_data_test_aux),2);

    % show the mean of the gaussian distribution
    GM.mu;
    GM.ComponentProportion;

    % take the mu with the maximim proportion
    %[~, index_max_mu] = max(GM.ComponentProportion);

    % order the phase in term of the proportion
    [~, index_sort] = sort(GM.ComponentProportion, "descend");

    mu_sorted = GM.mu(index_sort);
    % take the phase in the middle
    phase_middle = mean(mu_sorted);

    % sum to the outliers pi/2
    if (mu_sorted(2) < mu_sorted(1))
        index_out = angle(csi_data_test_aux) < phase_middle;
    else
        index_out = angle(csi_data_test_aux) > phase_middle;
    end

    % check whether is pi/2, pi and 3pi/2
    grid_pi = [pi/2, pi, -pi/2];

    diff_mu = diff(mu_sorted);

    % look for the minimum of the substraction
    sub_diff_mu = diff_mu + grid_pi;
    [~, index_min_mu] = min(abs(sub_diff_mu));


    % this is hard-coding --> analizy when the jump is pi/2, pi or 3*pi/2. this
    % case is 3pi/2 which is equal to -pi/2
    csi_data_test_aux(index_out) = csi_data_test_aux(index_out) .* exp(1i*grid_pi(index_min_mu));

    % csi_data_test_aux = csi_data_test .* exp(-1i*offset_mu);


    % remove the pi/2 offset
    csi_data_test_aux = csi_data_test_aux .* exp(-1i*pi/2);

    %figure, plot(angle(csi_data_test_aux));

    % it should be center at 0 or closer to 0
    mean(angle(csi_data_test_aux));

    % csi_data_test = csi_data_test*conj(exp(-1i*grid(index_max));
    % figure, plot(angle(csi_data_test));

    % figure, plot(grid,real(sum(spectrum_phase,1)))
    
    offset_phase = offset_phase - mean(angle(csi_data_test_aux));


end
