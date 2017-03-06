function [rb_idx, rb_mean, rb_ci, stat] = data_mse_analysis(P, s, d_snr, snr_hat, P_hat)
% DATA_MSE_ANALYSIS analysis of mean-square error E(\|s-\hat{s}\|^2)
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

% Copyright 2017 Richard J. Cui. Created: Mon 02/27/2017  2:15:49.215 PM
% $Revision: 0.3 $  $Date: Thu 03/02/2017 11:11:42.639 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% mean square error
% =========================================================================
[num_test, num_snr] = size(snr_hat);
rs = real(s);
N = length(rs);
Q = size(P, 1);

err_mpem = get_serr(num_test, num_snr, rs, N, Q, P_hat, 'MPEM');
err_mle  = get_serr(num_test, num_snr, rs, N, Q, P_hat, 'MLE');

% =========================================================================
% statistics
% =========================================================================
para.nboot = 1000;
para.nSNR = num_snr;
para.BaseIndex = find(d_snr == Inf);
[rb_idx, rb_mean, rb_ci] = get_robust_index(err_mpem, err_mle, para); % robustness index of MPEM algorithm

% ----------------------
% significance test
% ----------------------
stat = sig_test(rb_idx);

end % function data_mse_analysis

% =========================================================================
% subroutines
% =========================================================================
function stat = sig_test(idx)

[M, N] = size(idx);
h = zeros(1, N);
p = zeros(1, N);
for k = 1:N
    [p(k), h(k)] = ranksum(idx(:, k), zeros(M, 1), 'tail', 'right');
end % for
stat = struct('h', h, 'p', p);

end % function

function z = bootstat(x, y)

mx = mean(x);
my = mean(y);

z = (my-mx)./(my+mx);

end % function

function [rb_idx, rb_m, rb_ci] = get_robust_index(err_mpem, err_mle, para)

base_idx = para.BaseIndex;
nboot = para.nboot;
nSNR = para.nSNR;

x = 1:nSNR;
mpem = err_mpem(:, ~(x == base_idx));
mle = err_mle(:,  ~(x == base_idx));

rng(100) % for repeatness
rb_idx = bootstrp(nboot, @bootstat, mpem, mle); % get robutness index
rb_m = mean(rb_idx);
rb_ci = bootci(nboot, {@bootstat, mpem, mle}, 'type', 'bca'); % get confidence interval

end % function

function err = get_serr(num_test, num_snr, rs, N, Q, P_hat, alg_type)
% calculate square error
% 
% alg_type      - 'MPEM' or 'MLE'

err = zeros(num_test, num_snr);

for k = 1:num_test
    for j = 1:num_snr
        % cal. square error
        Phat_kj = zeros(Q, 5);
        Phat_kj(:, 1) = P_hat.(alg_type).A(k, j, :);
        Phat_kj(:, 2) = P_hat.(alg_type).tc(k, j, :);
        Phat_kj(:, 3) = P_hat.(alg_type).fc(k, j, :);
        Phat_kj(:, 4) = P_hat.(alg_type).cr(k, j, :);
        Phat_kj(:, 5) = P_hat.(alg_type).d(k, j, :);     
        
        rs_hat_kj = real(make_chirplets(N, Phat_kj));
        err(k, j) = norm(rs - rs_hat_kj)^2;        
    end % for
end % for


end % function

% [EOF]
