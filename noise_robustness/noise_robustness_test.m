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
d_snr       = [-1.0, 0.0]; % central values of desired SNR
num_test    = 2; % number of test at each test point

% =========================================================================
% test MPEM-ACT and MLE algorithm and get the data
% =========================================================================
Q = size(P, 1); % number of chiprlets to be estimated
[snr_hat, tcost, P_hat] = noise_robust_test(s, Q, d_snr, num_test);

% =========================================================================
% data analysis and visualization
% =========================================================================


% [EOF]
