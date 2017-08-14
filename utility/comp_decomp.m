function comp_decomp(s, spn, P, varargin)
% COMP_DECOMP compare original and recontructed signals from decomposition
%
% Syntax:
%   comp_decomp(s, spn, P)
%   comp_decomp(____, fig_name)
%   comp_decomp(____, 'PType', p_type)
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
% $Revision: 0.2 $  $Date: Wed 05/31/2017  4:35:55.161 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

p = check_inputs(s, spn, P, varargin{:});
s = p.s;
spn = p.spn;
P = p.P;
fig_name = p.fig_name;
p_type = p.PType;

N = length(spn);
rs = real(make_chirplets(N, P, 'PType', p_type)); % reconstructed signal

switch p_type
    case 'cohen'
        fs = N;
    case 'oneill'
        fs = 1;
    otherwise
        fs = 1;
end % switch
t = (0:N-1)/fs;

% compare original & reconstructed signals
% ----------------------------------------
figure('Name', fig_name)
subplot(311)
hs = plot(t, s);
axis tight, grid on;
hold on
hns = plot(t, spn);
hold off
legend([hs, hns], {'Clean' 'Noisy'})
ax_lim = axis(gca);
title 'Clean and noisy signal';

subplot(312)
hns = plot(t, spn);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(t, rs);
hold off
legend([hns, hrs], {'Noisy' 'Recon'})
title 'Noisy and reconstructed signal';

subplot(313), plot(t, spn-rs), axis(gca, ax_lim), grid on;
xlabel('Time (s)')
title 'Residuals = Noisy - Recon';

% compare clean & reconstructed signals
% ----------------------------------------
figure('Name', fig_name)
subplot(311)
hs = plot(t, s);
axis tight, grid on;
hold on
hns = plot(t, spn);
hold off
legend([hs, hns], {'Clean' 'Noisy'})
ax_lim = axis(gca);
title 'Clean and noisy signal';

subplot(312)
hs = plot(t, s);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(t, rs);
hold off
legend([hs, hrs], {'Clean' 'Recon'})
title 'Clean and reconstructed signal';

subplot(313), plot(t, s - rs), axis(gca, ax_lim), grid on;
xlabel('Time (s)')
title 'Error = Clean - Recon';

end % function comp_decomp

% =========================================================================
% subroutines
% =========================================================================
function q = check_inputs(s, spn, P, varargin)

ptype_str = {'oneill', 'cohen'};

p = inputParser;

p.addRequired('s', @isnumeric);
p.addRequired('spn', @isnumeric);
p.addRequired('P', @isnumeric);
p.addOptional('fig_name', '', @ischar);
p.addParameter('PType', 'Oneill', @(x) any(validatestring(lower(x), ptype_str)));

p.parse(s, spn, P, varargin{:});

q = p.Results;
q.PType = lower(q.PType);

end % function

% [EOF]
