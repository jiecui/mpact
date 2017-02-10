function varargout = mp_act_signal(s, Q, M, D, i0, a, ref_alg, verbose, mnits, level)
% MP_ACT_SIGNAL perform MP adaptive chirplet decomposition and show the results
%
% Syntax:
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

% Copyright 2016 Richard J. Cui. Created: Sat 12/17/2016 10:33:10.954 AM
% $Revision: 0.1 $  $Date: Sat 12/17/2016 10:33:10.965 AM $
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
show_decomp(s, P, fs, ref_alg)

% output
% ------
if nargout > 0
    varargout{1} = P;
end % if

end % function mp_act_signal

% [EOF]
