function [P, e, res] = em_chirplets(x, P, res, M, D, i0, radix, mnits,...
    verbose, mstep)
% EM_CHIRPLETS refine multiple chirplets with Expectation-Maximization algorithm
%
% Syntax:
%   [P, e, res] = em_chirplets(x, P, res, M, D, i0, radix, mnits, verbose)
% 
% Input(s):
%   x           - signal
%   P           - vector of chirplet parameters (see make_chirplets.m)
%   res         - residual = norm(signal - chirplets)
%   M           - resolution for Newton-Raphson refinement 
%   radix       - radix of scale
%   mnits       - maximum number of iterations
%   verbose     - verbose flag
%   mstep       - method of m-step, 'ONEILL' (default) or 'FEDER'
%
% Output(s):
%   P           - vector of chirplet parameters (see make_chirplets.m)
%   e           - e = signal - chirplets
%   res         - residual = norm(signal - chirplets)
%
% Example:
%
% Note:
%   See mp_adapt_chirplets.m for some explanation of the input/output
%   variables.
%
% References:
%   [1]	S. Mann and S. Haykin, "The adaptive chirplet - an adaptive 
%       generalized wavelet-like transform," in Adaptive Signal Processing.
%       vol. 1565, S. Haykin, Ed., ed Bellingham: SPIE - Int Soc Optical
%       Engineering, 1991, pp. 402-413.
%   [2]	J. Cui, "Adaptive chirplet transform for the analysis of visual 
%       evoked potentials," Doctor of Philosophy Dissertation (Ph.D.),
%       University of Toronto, 2006.
%   [3]	M. Feder and E. Weinstein, "Parameter estimation of superimposed
%       signals using the EM algorithm," IEEE Transactions on Acoustics,
%       Speech, and Signal Processing, vol. 36, pp. 477-489, 1988/04//
%       1988.
%   [4]	J. C. O'Neill and P. Flandrin, "Chirp hunting," in Proceedings of 
%       the IEEE-SP International Symposium on Time-Frequency and
%       Time-Scale Analysis, 1998, pp. 425-428.
% 
% See also make_chirplets, mp_adapt_chirplets.

% Copyright 2016-2017 Richard J. Cui. Created: Wed 12/14/2016  9:47:28.260 PM
% $Revision: 0.5 $  $Date: Sat 05/06/2017  2:33:57.838 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% Input parameters and options
% =========================================================================
if ~exist('mstep', 'var')
    mstep = 'oneill';
end % if

N   = length(x); % signal length
Q   = size(P, 1); % number of chirplets
P0  = zeros(Q, 5); % initial values of chirplet parameters
Ts  = ones(Q,1)*[.1 .001 1e-5 .1]; % tolerance of parameter error

% =========================================================================
% expectation-maximization (EM)
% =========================================================================
j = 1;
Pe = abs(P(:, 2:5) - P0(:, 2:5)); % parameter error
while sum(sum(Pe > Ts)) && (j <= mnits) && (Q > 1)
    P0 = P;
    
    % ------
    % E-Step
    % ------
    z = make_chirplets(N, P0);
    d = x - z; % note: d should be cloase to WGN
    
    % ------
    % M-Step
    % ------
    P = cal_mstep(j, N, P, d, M, D, i0, radix, mstep);
    Pe = abs(P(:, 2:5) - P0(:, 2:5));
    if verbose, fprintf('%d ', j); end
    j = j + 1;
end % while

% =========================================================================
% update the residual
% =========================================================================
y = make_chirplets(N,P);
e = x - y;
res = cat(1, res, norm(e));

end % function em_chirplets

% =========================================================================
% update the residual
% =========================================================================
function P = mstep_oneill(j, N, P, d, M, D, i0, radix)
% O'Neill method of m-step

Q = size(P, 1); % number of chirplets

for k = 1:Q
    if mod(k, j) == 0 % O'Neill's method
        z_k = make_chirplets(N, P(k,:));
        y_k = z_k + d;
        P_k = mp_chirplet(y_k, M, D, i0, radix, false);
        P(k,:) = P_k;
    end % if
end % for

end % function

function P = mstep_feder(N, P, d, M, D, i0, radix)
% Feder method of m-step

Q = size(P, 1); % number of chirplets

for k = 1:Q
    z_k = make_chirplets(N, P(k,:));
    y_k = z_k + d/Q;
    P_k = mp_chirplet(y_k, M, D, i0, radix, false);
    P(k,:) = P_k;
end % for

end % function

function P = cal_mstep(j, N, P, d, M, D, i0, radix, mstep)
% calculate m-step

switch lower(mstep)
    case 'feder' % refine all the chirplets at each iteration
        P = mstep_feder(N, P, d, M, D, i0, radix);
    case 'oneill' % refine only one chirplet at each iteration
        P = mstep_oneill(j, N, P, d, M, D, i0, radix);
end % switch

end % function

% [EOF]
