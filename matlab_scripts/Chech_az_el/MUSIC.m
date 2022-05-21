function [ps_db, D] = MUSIC(C, S, n_s, enableGPu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    % The number of samples
    
    if nargin < 3
        enableGPu = false;
        n_e = 1;
    end


    % Estimate eigen values and vectors by eigen decomposition
    % U contains the eigen vector and D contains the eigen value
    [U,D] = eig(C);
    % take the diagonal
    D = diag(D);
    % sort them by power
    [D,ind] = sort(D, 'descend');
    U = U(:,ind); 

    % It is needed to estimate the noise space, taking into account that we
    % have n_s signals
    Un = U(:,(n_s + 1):end);

    % try to do it by GPU if not do it by CPU
    try
    canUseGPU = parallel.gpu.GPUDevice.isAvailable;
    catch ME
    canUseGPU = false;
    end
    if(canUseGPU && enableGPu)
        S_gpu = gpuArray(S);
        Un_gpu = gpuArray(Un);
        ps = sum(abs((S_gpu')*Un_gpu).^2,2);
        ps = ps.^-1;
        ps = gather(ps);
    else
        % Pseudo-spectrum
%         tic
        ps = sum(abs((S')*Un).^2,2);
%         ps = real(ps);
        ps = ps.^-1;
%         toc
    end

    % conver to log scale
    ps_db = 10*log10(ps);

end

