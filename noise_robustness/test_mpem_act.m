function [P, t] = test_mpem_act(Q, s, varargin)
% TEST_MPEM_ACT test noise robustness of MPEM-ACT algorithm
%
% Syntax:
%   P = test_mpem_act(Q, spn)
%   P = test_mpem_act(____, 'PType', p_type)
%   P = test_mpem_act(____, 'Verbose', verb)
%   [P, t] = test_mpem_act(____)
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
% Example:
%
% Note:
%
% References:
%
% See also make_chirplets.

% Copyright 2017 Richard J. Cui. Created: Wed 02/22/2017  3:38:33.313 PM
% $Revision: 0.2 $  $Date: Wed 05/31/2017  1:18:33.652 PM $
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

p = do_mpem_act(Q, s, p_type, verb);
P = array2table(p, 'VariableNames', {'A', 'tc', 'fc', 'cr', 'd'});

t = toc(tv);

end % function test_mpem_act

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

function P = do_mpem_act(Q, spn, p_type, verbose)
% chirplet decomposition with MPEM algorithm

% required parameters
% --------------------
i0  = 1; % the first scale to roate the atoms
D   = 5; % decomposition depth = the higest scale
a   = 2; % the radix of scale
M = 256; % resolution for Newton-Raphson refinement
mnits   = 10; % max number of iteration for refinement
level   = 2; % diffculty level for MLE refinement

P = mp_adapt_chirplets(spn, Q, M, D, i0, a, verbose, mnits, level,...
    'RefineAlgorithm', 'ExpectMax', 'PType', p_type);

end % function

% [EOF]
