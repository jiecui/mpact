function [snr_hat, tcost, p_hat] = noise_robust_test(s, Q, d_snr, num_test)
% NOISE_ROBUST_TEST test noise robustness of MPEM-ACT and MLE
%
% Syntax:
%
% Input(s):
%   s           - clean (complex) signal
%   Q           - number of chirplets desired
%   d_snr       - 1 x N desired SNR
%   num_test    - number of tests for each d_snr
% 
% Output(s):
%   snr_hat     - estimated SNR for each test, num_dsnr x num_test array
%   tcost       - time cost (same structure as snr_hat)
%                 .MPEM
%                 .MLE
%   p_hat       - estimated chirplet parameters
%                 .MPEM (A, tc, fc, cr, d)
%                 .MLE  (A, tc, fc, cr, d)
%                 num_dsnr x num_test array x Q
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2017 Richard J. Cui. Created: Tue 02/21/2017  3:14:48.375 PM
% $Revision: 0.1 $  $Date: Tue 02/21/2017  3:14:48.380 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% main
% =========================================================================
num_dsnr = numel(d_snr);
x = zeros(num_test, num_dsnr); 
y = zeros(num_test, num_dsnr, Q);
z = struct('A', y, 'tc', y, 'fc', y, 'cr', y, 'd', y);

snr_hat = x;
tcost = struct('MPEM', x, 'MLE', x);
p_hat = struct('MPEM', z, 'MLE', z);

hw = waitbar(0, 'Noise robustness testing ...');
for m = 1:num_test
    waitbar((m-1)/num_test, hw)
    for n = 1:num_dsnr
        % get noisy signal & estimated SNR
        % --------------------------------
        dsnr_n = d_snr(n); % nth desired SNR
        [spn_mn, ~, snr_hat(m, n)] = add_noise(s, dsnr_n); % add noise to the signal
        
        % test MPEM-ACT algorithm
        % -----------------------
        [P_mn, t_mn] = test_mpem_act(Q, spn_mn);
        tcost.MPEM(m, n)        = t_mn;
        p_hat.MPEM.A(m, n, :)   = P_mn.A;
        p_hat.MPEM.tc(m, n, :)  = P_mn.tc;
        p_hat.MPEM.fc(m, n, :)  = P_mn.fc;
        p_hat.MPEM.cr(m, n, :)  = P_mn.cr;
        p_hat.MPEM.d(m, n, :)   = P_mn.d;
        
        % test MLE algorithm
        % ------------------
        [P_mn, t_mn] = test_mle_act(Q, spn_mn);
        tcost.MLE(m, n)         = t_mn;
        p_hat.MLE.A(m, n, :)    = P_mn.A;
        p_hat.MLE.tc(m, n, :)   = P_mn.tc;
        p_hat.MLE.fc(m, n, :)   = P_mn.fc;
        p_hat.MLE.cr(m, n, :)   = P_mn.cr;
        p_hat.MLE.d(m, n, :)    = P_mn.d;
    end % for
end % for
waitbar(m/num_test, hw), pause(.1)
close(hw)

end % function noise_robust_test

% [EOF]
