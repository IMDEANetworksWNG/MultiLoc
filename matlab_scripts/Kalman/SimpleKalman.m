function [X, sigma_x] = SimpleKalman(Y, sigma_m, sigma_e)

    X=Y;
    sigma_x=repmat(sigma_m, 1, size(Y, 2));
    sigma_m_inv=1/sigma_m;

    for ii=2:length(sigma_x)
        
        sigma_x(ii)=(sigma_x(ii-1)^2+sigma_e^2)^(-1/2);
        X(:, ii)=(sigma_x(ii)+sigma_m_inv)^-1*(sigma_x(ii)*X(:, ii-1)+sigma_m_inv*Y(:, ii));
        sigma_x(ii)=(sigma_x(ii)^2+sigma_m_inv^2)^(-1/2);
        
    end

end