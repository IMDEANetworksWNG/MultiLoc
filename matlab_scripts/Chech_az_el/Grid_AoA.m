function [cb_aoa, theta] = Grid_AoA(step_angle, N,d,lambda)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    
    %%%%%% AoA %%%%%%
    theta = -90:step_angle:(90-step_angle);
    % pass to radian
    theta = deg2rad(theta);

    % Calculate the steering matrix
    cb_aoa = ones(N,length(theta));

    for i = 1:N
        cb_aoa(i,:) = exp(-1i*(i-1)*d*2*pi*sin(theta)/lambda);
    end
end

