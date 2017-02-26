% NOISE_ROBUSTNESS_TEST Test noise robustness of MPEM and MLE algorithms

% Copyright 2017 Richard J. Cui. Created: Mon 02/20/2017  2:14:47.141 PM
% $Revision: 0.1 $  $Date: Mon 02/20/2017  2:14:47.149 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% set parameters
% =========================================================================
% chirplet signals
% ----------------
N = 100; % signal length
P = [1*exp(1i*0), N/2+1, pi/2,  pi/N, N/3; % up-chirplet
     1*exp(1i*0), N/2+1, pi/2, -pi/N, N/3]; % down-chirplet
s = make_chirplets(N, P); % the complex signal

% test paras
% ----------
d_snr       = [-3.0, Inf]; % central values of desired SNR
num_test    = 2; % number of test at each test point

% =========================================================================
% test MPEM-ACT and MLE algorithm and get the data
% =========================================================================
Q = size(P, 1); % number of chiprlets to be estimated
[snr_hat, tcost, P_hat] = noise_robust_test(s, Q, d_snr, num_test);

% =========================================================================
% data analysis and visualization
% =========================================================================
%% ---------------------------------------
% mean square error = E(\|s-\hat{s}\|^2)
% ---------------------------------------
rs = real(s);
num_snr = numel(d_snr);
err_mpem = zeros(num_test, num_snr);
err_mle = err_mpem;
for k = 1:num_test
    for j = 1:num_snr
        % MPEM
        PhatMpem_kj = zeros(Q, 5);
        PhatMpem_kj(:, 1) = P_hat.MPEM.A(k, j, :);
        PhatMpem_kj(:, 2) = P_hat.MPEM.tc(k, j, :);
        PhatMpem_kj(:, 3) = P_hat.MPEM.fc(k, j, :);
        PhatMpem_kj(:, 4) = P_hat.MPEM.cr(k, j, :);
        PhatMpem_kj(:, 5) = P_hat.MPEM.d(k, j, :);     
        
        rs_hat_kj = real(make_chirplets(N, PhatMpem_kj));
        err_mpem(k, j) = norm(rs - rs_hat_kj)^2;
        
        % MLE
        PhatMle_kj = zeros(Q, 5);
        PhatMle_kj(:, 1) = P_hat.MLE.A(k, j, :);
        PhatMle_kj(:, 2) = P_hat.MLE.tc(k, j, :);
        PhatMle_kj(:, 3) = P_hat.MLE.fc(k, j, :);
        PhatMle_kj(:, 4) = P_hat.MLE.cr(k, j, :);
        PhatMle_kj(:, 5) = P_hat.MLE.d(k, j, :);     
        
        rs_hat_kj = real(make_chirplets(N, PhatMle_kj));
        err_mle(k, j) = norm(rs - rs_hat_kj)^2;        
    end % for
end % for


% [EOF]
