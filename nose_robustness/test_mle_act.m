function [t, P] = test_mle_act(Q, spn)
% TEST_MLE_ACT test noise robustness of MLE-ACT algorithm
%
% Syntax:
%
% Input(s):
%   Q           - number of chiprlets desired
%   spn         - signal + noise
%
% Output(s):
%   t           - time cost of this run (seconds)
%   P           - [A, tc, fc, cr, d; ...]
%                   A   :amplitude of the estimated chirplets 1 x Q array
%                   tc  :time-center
%                   fc  :frequency-center
%                   cr  :chirp-rate
%                   d   :duration
%
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

% Copyright 2017 Richard J. Cui. Created: Thu 02/23/2017  3:31:37.681 PM
% $Revision: 0.1 $  $Date: Thu 02/23/2017  3:31:37.694 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% main
% =========================================================================
tv = tic;

p = do_mle_act(Q, spn);
P = array2table(p, 'VariableNames', {'A', 'tc', 'fc', 'cr', 'd'});

t = toc(tv);

end % function test_mle_act

% =========================================================================
% subroutines
% =========================================================================
function P = do_mle_act(Q, spn)
% chirplet decomposition with MLE algorithm

% Common parameters
% -----------------
d   = 50; % initial guess of duration
cr  = 0; % chirprate guess
M   = 256; % resolution for Newton-Raphson refinement
verbose = 'No'; % don't show notes
mnits   = 10; % max number of iteration for refinement
level   = 2; % difficult level

P = mle_adapt_chirplets(spn, Q, M, d, cr, verbose, mnits, level,...
    'RefineAlgorithm', 'ExpectMax');

end % function

% [EOF]
