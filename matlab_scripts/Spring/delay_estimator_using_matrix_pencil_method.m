function [all_gains, all_delays, num_estimated_paths,matrix_estimated_gains,matrix_estimated_delay] = delay_estimator_using_matrix_pencil_method(multipath_channel,P,z_left,z_right,unit_circle_toll)
    [antennas, N] = size(multipath_channel);
    num_subcarriers = N-z_left-z_right;
    all_gains = NaN.*ones(antennas, 100);
    all_delays = NaN.*ones(antennas, 100);
    num_estimated_paths = zeros(1,antennas);
    matrix_estimated_gains = [];
    matrix_estimated_delay = [];
    for antenna = 1:antennas
        H_n = transpose(multipath_channel(antenna,z_left+1:end-z_right));
        % Hankel matrix:
        K = num_subcarriers-P;
        Hankel_matrix = [];
        for k = 1:K
            Hankel_matrix(k,:) = H_n(k:k+P);
        end
        H1 = Hankel_matrix(:,1:P);
        H2 = Hankel_matrix(:,2:P+1);
        H_H = pinv(H2)*H1;
        z_l = eig(H_H);
        
        estimated_delay = (N/(2*pi)).*(angle(z_l));
        
        %----distance to the unit circle----% 
        d = abs(1-abs(z_l));
        %-----------------------------------%
        
        ind_d = find(d<unit_circle_toll);
        if isempty(ind_d)
            [~, ind_d] = min(d); 
        end
        delays = sort(estimated_delay(ind_d));
        A_tau = [];    
        for n = 1:num_subcarriers
            z_l_per_subcarrier = zeros(1, length(delays));
           for l = 1:length(delays)
              z_l_per_subcarrier(l) = exp(-1i*2*pi*delays(l)*n/N); 
           end
           A_tau = vertcat(A_tau, z_l_per_subcarrier);
        end
        gains = pinv(A_tau) * H_n;
        all_gains(antenna,1:length(gains)) = transpose(gains);
        all_delays(antenna,1:length(gains)) = transpose(delays);
        num_estimated_paths(antenna) = length(gains);
    end
end

