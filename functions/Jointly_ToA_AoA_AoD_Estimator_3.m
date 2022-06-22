function [matrix_toa_aoa_aod] = Jointly_ToA_AoA_AoD_Estimator_3(S_toa,S_aoa,S_aod, channel)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[~,len_aod] = size(S_aod);
[~,len_aoa] = size(S_aoa);
[~,len_toa] = size(S_toa);
matrix_toa_aoa_aod = zeros(len_toa, len_aoa, len_aod);
%     matrix_toa_aoa_aod = zeros(len_toa, len_aoa);

[K,N,M] = size(channel);

% look over the aod
channel_reshape = reshape(channel, K*N, M);
channel_aod = channel_reshape * conj(S_aod);
channel_aod = reshape(channel_aod, K, N, len_aod);

% look over the aoa
channel_aoa = permute(channel_aod, [1,3,2]);
channel_aoa = reshape(channel_aoa, K*len_aod, N);
channel_aoa = channel_aoa * conj(S_aoa);
channel_aoa = reshape(channel_aoa, K,len_aod, len_aoa);
channel_aoa = permute(channel_aoa, [1,3,2]);

% look over the toa
channel_toa = permute(channel_aoa, [2,3,1]);
channel_toa = reshape(channel_toa, len_aoa*len_aod, K);
channel_toa = channel_toa * conj(S_toa);
channel_toa = reshape(channel_toa, len_aoa, len_aod, len_toa);

matrix_toa_aoa_aod = permute(channel_toa, [3,1,2]);
matrix_toa_aoa_aod = matrix_toa_aoa_aod./(K*N*M);
%     channel_aod = sum(channel_aod,3)/M;
%     for index_aoa = 1:len_aoa
%         channel_aoa = channel_aod .* S_aoa(:,index_aoa)';
%         channel_aoa = sum(channel_aoa,2)/N;
% 
%         channel_toa = (S_toa')*channel_aoa;
%         alpha = channel_toa/K;
% 
% %         power = alpha.^2;
%         matrix_toa_aoa_aod(:,index_aoa, index_aod) = alpha;
%     end
% end
end

