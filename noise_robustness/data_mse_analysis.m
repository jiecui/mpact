function [e_mpem, e_mle] = data_mse_analysis(P, s, d_snr, snr_hat, P_hat)
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
% $Revision: 0.1 $  $Date: Tue 02/28/2017 11:41:58.340 AM $
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
para.BaseIndex = find(d_snr == Inf);
e_mpem = err_m_sem(err_mpem, para); % find mean $ SEM of square-error
e_mle = err_m_sem(err_mle, para);

end % function data_mse_analysis

% =========================================================================
% subroutines
% =========================================================================
function m_sem = err_m_sem(serr, para)
% mean and sem of error

nboot = para.nboot;
rbi = get_robust_index(serr, para); % robustness index of MPEM algorithm
m_sem = [mean(rbi)', std(rbi)'/sqrt(nboot)]; % [mean, sem]

end % function

function rb_idx = get_robust_index(err_mse, para)

base_idx = para.BaseIndex;
nboot = para.nboot;

x = 1:size(err_mse, 2);
E = err_mse(:, ~(x == base_idx));
% B = err_mse(:, base_idx)*ones(1, size(E, 2));

err_idx = E;
% err_idx = (E-B)./(E+B);
% err_idx = (E-B)./B;

rb_idx = bootstrp(nboot, @(x) mean(x), err_idx); % get robutness index

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
