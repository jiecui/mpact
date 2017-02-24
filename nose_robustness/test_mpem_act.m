function [t, P] = test_mpem_act(Q, spn)
% TEST_MPEM_ACT test noise robustness of MPEM-ACT algorithm
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
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2017 Richard J. Cui. Created: Wed 02/22/2017  3:38:33.313 PM
% $Revision: 0.1 $  $Date: Wed 02/22/2017  3:38:33.317 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% main
% =========================================================================
tv = tic;

p = do_mpem_act(Q, spn);
P = array2table(p, 'VariableNames', {'A', 'tc', 'fc', 'cr', 'd'});

t = toc(tv);

end % function test_mpem_act

% =========================================================================
% subroutines
% =========================================================================
function P = do_mpem_act(Q, spn)
% chirplet decomposition with MPEM algorithm

% Common parameters
% -----------------
i0  = 1; % the first scale to roate the atoms
D   = 5; % decomposition depth = the higest scale
a   = 2; % the radix of scale
M = 256; % resolution for Newton-Raphson refinement
verbose = 'No'; % don't show notes
mnits   = 10; % max number of iteration for refinement
level   = 2; % diffculty level for MLE refinement

P = mp_adapt_chirplets(spn, Q, M, D, i0, a, verbose, mnits, level,...
    'RefineAlgorithm', 'ExpectMax');

end % function

% [EOF]
