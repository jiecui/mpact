function comp_sqerr(s, P_mpem, P_mle)
% COMP_SQERR compare square error of MPEM and MLE
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

% Copyright 2017 Richard J. Cui. Created: Wed 03/08/2017 10:32:27.292 AM
% $Revision: 0.2 $  $Date: Wed 03/08/2017 10:32:27.292 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% Input parameters and options
% =========================================================================
fig_name = 'MPEM vs. MLE';
N = length(s);
p_mpem = table2array(P_mpem);
p_mle = table2array(P_mle);
rs_mpem = real(make_chirplets(N, p_mpem)); % reconstructed mpem signal
rs_mle  = real(make_chirplets(N, p_mle)); % reconstructed mle signal

err_mpem = abs(s-rs_mpem).^2;
err_mle  = abs(s-rs_mle).^2;

% =========================================================================
% compare original & reconstructed signals
% =========================================================================
figure('Name', fig_name)
subplot(311) % original and recon with MPEM
hs = plot(s);
axis tight, grid on;
hold on
hrs = plot(rs_mpem);
hold off
legend([hs, hrs], {'Clean' 'Recon'})
ax_lim = axis(gca);
title 'Clean and MPEM Reconstructed signal';

subplot(312) % original and recon with MLE
hs = plot(s);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(rs_mle);
hold off
legend([hs, hrs], {'Clean' 'Recon'})
title 'Clean and MLE reconstructed signal';

subplot(313) % compare square error
hmpem = plot(err_mpem);
hold on
hmle = plot(err_mle);
hold off
grid on
legend([hmpem, hmle], {'MPEM', 'MLE'})
xlabel('Time (s)')
title 'Compare squared error';

end % function comp_decomp

% [EOF]
