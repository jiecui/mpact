function comp_sqerr(s, p_mpem, p_mle, varargin)
% COMP_SQERR compare square error of MPEM and MLE
%
% Syntax:
%   comp_sqerr(s, P_mpem, P_mle)
%   comp_sqerr(____, 'PType', p_type)
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
p = check_inputs(s, p_mpem, p_mle, varargin{:});
p_mpem = p.p_mpem;
p_mle = p.p_mle;
p_type = p.PType;

N = length(s);
fig_name = 'MPEM vs. MLE';
switch p_type
    case 'cohen'
        fs = N;
    case 'oneill'
        fs = 1;
    otherwise
        fs = 1;
end % switch
t = (0:N-1)/fs;

rs_mpem = real(make_chirplets(N, p_mpem, 'PType', p_type)); % reconstructed mpem signal
rs_mle  = real(make_chirplets(N, p_mle, 'PType', p_type)); % reconstructed mle signal

err_mpem = abs(s-rs_mpem).^2;
err_mle  = abs(s-rs_mle).^2;

% =========================================================================
% compare original & reconstructed signals
% =========================================================================
figure('Name', fig_name)
subplot(311) % original and recon with MPEM
hs = plot(t, s);
axis tight, grid on;
hold on
hrs = plot(t, rs_mpem);
hold off
legend([hs, hrs], {'Clean' 'Recon'})
ax_lim = axis(gca);
title 'Clean and MPEM Reconstructed signal';

subplot(312) % original and recon with MLE
hs = plot(t, s);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(t, rs_mle);
hold off
legend([hs, hrs], {'Clean' 'Recon'})
title 'Clean and MLE reconstructed signal';

subplot(313) % compare square error
hmpem = plot(t, err_mpem);
hold on
hmle = plot(t, err_mle);
hold off
grid on
legend([hmpem, hmle], {'MPEM', 'MLE'})
xlabel('Time (s)')
title 'Compare squared error';

end % function comp_decomp

% =========================================================================
% subroutines
% =========================================================================
function q = check_inputs(s, p_mpem, p_mle, varargin)

ptype_str = {'oneill', 'cohen'};

p = inputParser;

p.addRequired('s', @isnumeric);
p.addRequired('p_mpem', @isnumeric);
p.addRequired('p_mle', @isnumeric);
p.addParameter('PType', 'Oneill', @(x) any(validatestring(lower(x), ptype_str)));

p.parse(s, p_mpem, p_mle, varargin{:});

q = p.Results;
q.PType = lower(q.PType);

end % function

% [EOF]
