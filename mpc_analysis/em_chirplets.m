function [P, e, res] = em_chirplets(x, P, res, M, D, i0, radix, mnits, verbose)
% EM_CHIRPLETS refine multiple chirplets with Expectation-Maximization algorithm
%
% Syntax:
%   [P, e, res] = em_chirplets(x, P, res, M, D, i0, radix, mnits, verbose)
% 
% Input(s):
%   x           - signal
%   P           - vector of chirplet parameters (see make_chirplets.m)
%   res         - residual = norm(signal - chirplets)
%   M           - resolution for Newton-Raphson refinement (optional, 
%   radix       - radix of scale (default = 2)
%   mnits       - maximum number of iterations (optional, default = 5)
%   verbose     - verbose flag (default = false)
%
% Output(s):
%   P           - vector of chirplet parameters (see make_chirplets.m)
%   e           - e = signal - chirplets
%   res         - residual = norm(signal - chirplets)
%
% Example:
%
% Note:
%
% References:
%   [1]	S. Mann and S. Haykin, "The adaptive chirplet - an adaptive 
%       generalized wavelet-like transform," in Adaptive Signal Processing.
%       vol. 1565, S. Haykin, Ed., ed Bellingham: SPIE - Int Soc Optical
%       Engineering, 1991, pp. 402-413.
%   [2]	J. Cui, "Adaptive chirplet transform for the analysis of visual 
%       evoked potentials," Doctor of Philosophy Dissertation (Ph.D.),
%       University of Toronto, 2006.
% 
% See also make_chirplets.

% Copyright 2016 Richard J. Cui. Created: Wed 12/14/2016  9:47:28.260 PM
% $Revision: 0.2 $  $Date: Fri 12/16/2016  3:11:09.454 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% Input parameters and options
% =========================================================================
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
    d = x - z;
    
    % ------
    % M-Step
    % ------
    for k = 1:Q
        z_k = make_chirplets(N, P(k,:));
        y_k = z_k + d; % Note: d/Q is not efficient; should look into it
        P_k = mp_chirplet(y_k, M, D, i0, radix, false);
        P(k,:) = P_k;
    end % for
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

end % function lem_chirplets

% [EOF]
