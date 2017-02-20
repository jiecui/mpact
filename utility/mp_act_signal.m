function varargout = mp_act_signal(s, Q, M, D, i0, a, ref_alg, verbose, mnits, level)
% MP_ACT_SIGNAL perform MP adaptive chirplet decomposition and show the results
%
% Syntax:
%
% Input(s):
%   ref_alg         - refinement algorithm
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

% Copyright 2016-2017 Richard J. Cui. Created: Sat 12/17/2016 10:33:10.954 AM
% $Revision: 0.3 $  $Date: Sun 02/19/2017 10:05:56.162 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

if ~exist('level', 'var')
    level = 2;
end % if

fs = 1; % assume normalized sampling frequency

% decompose signal
% ----------------
x = hilbert(s); % convert the signal into an analytical signal
P = mp_adapt_chirplets(x, Q, M, D, i0, a, verbose, mnits, level,...
    'RefineAlgorithm', ref_alg);

% Display the results of decomposition
% -----------------------------------
% get figure name
fig_name = get_fig_name(ref_alg);
% show figure
show_decomp(s, P, fs, fig_name)

% output
% ------
if nargout > 0
    varargout{1} = P;
end % if

end % function mp_act_signal

% [EOF]
