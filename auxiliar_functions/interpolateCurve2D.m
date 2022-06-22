function [X_interpolated, Y_interpolated] = interpolateCurve2D(X, Y, n)

    traj_x = X;
    traj_y = Y;

    traj_y_1 = traj_y;
    traj_y_1 = traj_y_1(1:(end-1));
    traj_y_2 = traj_y;
    traj_y_2 = traj_y_2(2:end);
    traj_y_interp = linspaceNDim(traj_y_1, traj_y_2, 100);
    traj_y_interp(:,end) = [];
    traj_y_interp = reshape(traj_y_interp.', numel(traj_y_interp),1);
    traj_y_interp(end+1) = traj_y(end);

    traj_x_1 = traj_x;
    traj_x_1 = traj_x_1(1:(end-1));
    traj_x_2 = traj_x;
    traj_x_2 = traj_x_2(2:end);
    traj_x_interp = linspaceNDim(traj_x_1, traj_x_2, 100);
    traj_x_interp(:,end) = [];
    traj_x_interp = reshape(traj_x_interp.', numel(traj_x_interp),1);
    traj_x_interp(end+1) = traj_x(end);
    
    % Now we have a bunch of points, pick the correct number of points
    bunch = size(traj_x_interp, 1);
    take_point_every_n_points = round(bunch/n);
    
    X_interpolated = traj_x_interp(1:take_point_every_n_points:end);
    Y_interpolated = traj_y_interp(1:take_point_every_n_points:end);

end

