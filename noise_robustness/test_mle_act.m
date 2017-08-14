function [P, t] = test_mle_act(Q, s, varargin)
% TEST_MLE_ACT test noise robustness of MLE-ACT algorithm
%
% Syntax:
%   P = test_mle_act(Q, spn)
%   P = test_mle_act(____, 'PType', p_type)
%   P = test_mle_act(____, 'Verbose', verb)
%   [P, t] = test_le_act(____)
% 
% Input(s):
%   Q           - number of chiprlets desired
%   s           - signal
%   PType       - (parameter) type of chirplet parameter P, p_type = 
%                 {'Oneill', 'Cohen'}
%   Verbose     - (Parameter) show details, verb = {'No', 
%                 'Yes', 'vv'} (default verb = 'No')
%
% Output(s):
%   P           - table of chirplet parameters [A, tc, fc, cr, d; ...]
%                   A   :amplitude of the estimated chirplets 1 x Q array
%                   tc  :time-center
%                   fc  :frequency-center
%                   cr  :chirp-rate
%                   d   :duration
%                 (see make_chirplets.m for the details)
%   t           - time cost of this run (seconds)
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
% $Revision: 0.2 $  $Date: Wed 05/31/2017  4:04:56.096 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% main
% =========================================================================
tv = tic;

p = check_inputs(Q, s, varargin{:});
Q = p.Q;
s = p.s;
p_type = p.PType;
verb = p.Verbose;

p = do_mle_act(Q, s,  p_type, verb);
P = array2table(p, 'VariableNames', {'A', 'tc', 'fc', 'cr', 'd'});

t = toc(tv);

end % function test_mle_act

% =========================================================================
% subroutines
% =========================================================================
function q = check_inputs(Q, s, varargin)

ptype_str = {'oneill', 'cohen'};
verb_str = {'No', 'Yes', 'vv'};

p = inputParser;

p.addRequired('Q', @isnumeric);
p.addRequired('s', @isnumeric);
p.addParameter('PType', 'Oneill', @(x) any(validatestring(lower(x), ptype_str)));
p.addParameter('Verbose', 'No', @(x) any(validatestring(lower(x), verb_str)));

p.parse(Q, s, varargin{:});

q = p.Results;
q.PType = lower(q.PType);
q.Verbose = lower(q.Verbose);

end % function

function P = do_mle_act(Q, spn, p_type, verbose)
% chirplet decomposition with MLE algorithm

% initial parameters
% -------------------
N   = length(spn); % signal length
tc  = N/2; % time-center
fc  = pi/2; % frequency-center
cr  = 0; % chirprate guess
d   = N/4; % initial guess of duration
M   = 256; % resolution for Newton-Raphson refinement
mnits   = 10; % max number of iteration for refinement
level   = 2; % difficult level

CP0 = [tc, fc, cr, d]; % initial chirplet parameters (no amplitude)
P = mle_adapt_chirplets(spn, Q, M, CP0, verbose, mnits, level,...
    'RefineAlgorithm', 'MaxLikeliEst', 'PType', p_type);

end % function

% [EOF]
