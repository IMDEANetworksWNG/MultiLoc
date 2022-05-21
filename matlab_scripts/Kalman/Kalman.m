function [x_kalman,y_kalman] = Kalman(x_noisy,y_noisy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %%%%% Kalman filter %%%%%

    power_noise = 1;
    % velocity. half meter per second
    v = 1.5;

    % initialize the target state
    state_init = [x_noisy(1), y_noisy(1), v, v];

    % allocate memory for updated states
    state_up = zeros(length(state_init),length(x_noisy));
    state_up(:,1) = state_init;

    % initilize the error covariance
    P = eye(length(state_init));


    % state every second
    delta_t = 1;

    % state transmission matrix
    STM =  [1 0 delta_t 0; ...
            0 1 0 delta_t; ...
            0 0 1     0  ; ...
            0 0 0      1  ];

    % error state transmission matrix
    %ro = [(delta_t^2)/2, (delta_t^2)/2, 1, 1]'.*eye(length(state));

    % noise state covarianza
    var_noise_state = 0.1;

    Q = eye(length(state_init))*var_noise_state;

    H = [1 0 0 0; 0 1 0 0];

    % covariance error observation
    R = eye(2)*power_noise;

    for ii = 1:(length(x_noisy)-1)

        %%%%%%%%%%% prediction %%%%%%%%%%%
        % state prediction
        state_pred = STM*state_up(:,ii);

        % state covariance prediction
        P_pred = STM*P*STM.' + Q;

        %%%%%%%%%%% update %%%%%%%%%%%
        % kalman gain
        K = P*H.'*inv(H*P*(H.') + R);

        % observation
        z = [x_noisy(ii);y_noisy(ii)];
        % state update
        state_up(:,ii+1) = state_pred + K*(z-H*state_pred);

        % state covariance update
        P = (eye(length(state_pred)) - K*H)*P_pred;
    end
    
    % take the coordinates
    x_kalman = state_up(1,:).';
    y_kalman = state_up(2,:).';
end

