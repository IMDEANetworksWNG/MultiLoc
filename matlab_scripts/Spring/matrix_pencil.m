function [z_l] = matrix_pencil(csi_data,pencil)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[sensors, aux] = size(csi_data);

if (aux > sensors)
    csi_data = csi_data.';
    [sensors, ~] = length(csi_data);

end


% Hankel matrix:
sub_k = sensors-pencil;

Hankel_matrix = [];
for k = 1:sub_k
    Hankel_matrix(k,:) = csi_data(k:k+pencil);
end

H = Hankel_matrix(:,1:(pencil+1));

% do the svd
[U,S,V] = svd(H);
S = diag(S);
V = conj(V);

% th for splitting the signal space
th = 1e-1;

% taje the max eigen-value
s_max = max(S);
% calculate the value of the min eigen value
s_min = s_max .* th;
% take the index
index_signal = S > s_min;

% take the first eigen right vectors
V_signal = V(:,index_signal);

V1 = V_signal(1:pencil,:);
V2 = V_signal(2:(pencil+1),:);
H_H = pinv(V2)*V1;

% take the eighen values
z_l = eig(H_H);
end

