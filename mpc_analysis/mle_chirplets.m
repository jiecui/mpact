function [P, e, res] = mle_chirplets(x, P, res, level, M, mnits, verbose)
% MLE_CHIRPLETS estimate chirplets with Maximum Likelihood approach
%
% Syntax:
%   [P, e, res] = mle_chirplets(x, P, res, M, emits, verbose)
% 
% Input(s):
%   x           - signal
%   P           - vector of chirplet parameters (see make_chirplets.m)
%   res         - residual = norm(signal - chirplets)
%   level       - level of difficulty of MLE
%   M           - resolution for Newton-Raphson refinement (optional, 
%   mnits       - maximum number of iterations (optional, default = 5)
%   verbose     - verbose flag (optional, default = false)
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
%  References
%   [1]	J. C. O'Neill and P. Flandrin, "Chirp hunting," in Proceedings of 
%       the IEEE-SP International Symposium on Time-Frequency and
%       Time-Scale Analysis, 1998, pp. 425-428.
% 
% See also .

% Copyright 2016 Richard J. Cui. Created: Wed 12/14/2016  9:47:28.260 PM
% $Revision: 0.1 $  $Date: Wed 12/14/2016  9:47:28.340 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% Input parameters and options
% =========================================================================
N = length(x);
Q   = size(P, 1); % number of chirplets
P0  = zeros(Q, 5); % initial values of chirplet parameters
Ts  = ones(Q,1)*[.1 .001 1e-5 .1]; % tolerance of parameter error

% =========================================================================
% Maximum-Likelihood Estimation (MLE)
% =========================================================================
j = 1;
Pe = abs(P(:, 2:5) - P0(:, 2:5)); % parameter error
while sum(sum(Pe > Ts)) && (j <= mnits) && (Q > 1)    
    % save current as previous
    P0 = P;
    z  = make_chirplets(N, P0);
    
    % update current with delta residual
    for k = 1:Q
        z_k = make_chirplets(N, P0(k,:));
        delta_k = x - (z - z_k);
        t_k = P(k, 2);
        f_k = P(k, 3);
        c_k = P(k, 4);
        d_k = P(k, 5);
        P_k = best_chirplet(delta_k, level, M, t_k, f_k, c_k, d_k, false);
        P(k,:) = P_k;
    end % for
    Pe = abs(P(:, 2:5) - P0(:, 2:5));
    if verbose, fprintf('%d ', j); end
    j = j + 1;
end % while

% get the residual
y = make_chirplets(N,P);
e = x - y;
res = cat(1, res, norm(e));
    
end % function lem_chirplets

% [EOF]
