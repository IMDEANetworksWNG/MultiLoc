function [matrix_toa_aoa_aod] = Jointly_ToA_AoA_AoD_Estimator_2(S_toa,S_aoa,S_aod, channel)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[~,len_aod] = size(S_aod);
[~,len_aoa] = size(S_aoa);
[~,len_toa] = size(S_toa);
matrix_toa_aoa_aod = zeros(len_toa, len_aoa, len_aod);
%     matrix_toa_aoa_aod = zeros(len_toa, len_aoa);

[K,N,M] = size(channel);

for index_aod = 1:len_aod
    S_aod_aux = reshape(S_aod(:,index_aod)',1,1,M);
    channel_aod = channel.* S_aod_aux;
    channel_aod = sum(channel_aod,3)/M;
    for index_aoa = 1:len_aoa
        channel_aoa = channel_aod .* S_aoa(:,index_aoa)';
        channel_aoa = sum(channel_aoa,2)/N;

        channel_toa = (S_toa')*channel_aoa;
        alpha = channel_toa/K;

%         power = alpha.^2;
        matrix_toa_aoa_aod(:,index_aoa, index_aod) = alpha;
    end
end
end

