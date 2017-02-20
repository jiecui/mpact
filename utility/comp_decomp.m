function comp_decomp(s, spn, P, fig_name)
% COMP_DECOMP compare original and recontructed signals from decomposition
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

% Copyright 2017 Richard J. Cui. Created: Mon 02/20/2017  9:32:38.830 AM
% $Revision: 0.1 $  $Date: Mon 02/20/2017  9:32:38.835 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

N = length(spn);
rs = real(make_chirplets(N, P)); % reconstructed signal

% compare original & reconstructed signals
% ----------------------------------------
figure('Name', fig_name)
subplot(311)
hs = plot(s);
axis tight, grid on;
hold on
hns = plot(spn);
hold off
legend([hs, hns], {'Clean' 'Noisy'})
ax_lim = axis(gca);
title 'Clean and noisy signal';

subplot(312)
hns = plot(spn);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(rs);
hold off
legend([hns, hrs], {'Noisy' 'Recon'})
title 'Noisy and reconstructed signal';

subplot(313), plot(spn-rs), axis(gca, ax_lim), grid on;
xlabel('Time (s)')
title 'Residuals = Noisy - Recon';

% compare clean & reconstructed signals
% ----------------------------------------
figure('Name', fig_name)
subplot(311)
hs = plot(s);
axis tight, grid on;
hold on
hns = plot(spn);
hold off
legend([hs, hns], {'Clean' 'Noisy'})
ax_lim = axis(gca);
title 'Clean and noisy signal';

subplot(312)
hs = plot(s);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(rs);
hold off
legend([hs, hrs], {'Clean' 'Recon'})
title 'Clean and reconstructed signal';

subplot(313), plot(s - rs), axis(gca, ax_lim), grid on;
xlabel('Time (s)')
title 'Error = Clean - Recon';

end % function comp_decomp

% [EOF]
