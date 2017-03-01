% SHOW_ROBUSTNESS_ANALYSIS display results
%
% Syntax:
%
% Input(s):
%
% Output(s):
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2017 Richard J. Cui. Created: Wed 03/01/2017  4:14:38.437 PM
% $Revision: 0.1 $  $Date: Wed 03/01/2017  4:14:38.441 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% data analysis
% =========================================================================
load noise_robust_data.mat
[e_mpem, e_mle] = data_mse_analysis(P, s, d_snr, snr_hat, P_hat);

% =========================================================================
% robustness indexes
% =========================================================================
f = @(s, e) errorbar(s(~(s == Inf)), e.Mean, e.Mean-e.CI(1, :), e.CI(2, :)-e.Mean);
figure
f(d_snr, e_mpem)
hold on
f(d_snr, e_mle)
ax = axis(gca);
xlim([ax(1)-5, ax(2)+5])
set(gca, 'YScale', 'log')
legend('MPEM', 'MLE')

% [EOF]
