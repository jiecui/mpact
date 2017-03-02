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
[rb_idx, rb_mean, rb_ci, stat] = data_mse_analysis(P, s, d_snr, snr_hat, P_hat);

% =========================================================================
% robustness indexes
% =========================================================================
f = @(s, m, c) errorbar(s(~(s == Inf)), m, m-c(1, :), c(2, :)-m, 'o');

x = d_snr(~(d_snr == Inf));
v = rb_mean;
xq = min(x):1:max(x);
vq = interp1(x, v, xq, 'pchip'); % Shape-preserving piecewise cubic interpolation

figure
eb = f(d_snr, rb_mean, rb_ci); % return errorbar object
hold on
plot(xq, vq)
ax = axis(gca);
xlim([ax(1)-5, ax(2)+5])
ylim([-.1 .4])
plot(xlim, [0 0])
legend(eb, 'Relative Rb')
xlabel('Tested SNR (dB)')
ylabel('Robustness Index')

% =========================================================================
% display significant test
% =========================================================================
cprintf('Keywords', 'Significant test at each SNR point:\n')
disp(stat)

% [EOF]
