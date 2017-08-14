function wig = chrpltwvd_direct(P, n, f)
% CHRPLTWVD_DIRECT compute chirplet spectrogram by estimating chriplet WVD directly
%
% Syntax:
%   wig = chrpltwvd_direct(P, n, f)
% 
% Input(s):
%   P       - matrix of parameters [amp time freq chirp_rate duration; ...]
%   n       - vector of time range (assume fs = 1 Hz)
%   f       - vector of frequency range ([0 1] range)
%
% Outputs:
%   wig     - length(n)_by_length(f) WVD
%
% Note:
%   Assume O'Neill's equation and the sampling frequency of chirplet is
%   normalized to 1 Hz.
% 
% Example:
%
% References:
%
% See also .

% Copyright 2016 Richard J. Cui. Created: Tue 12/13/2016 11:58:46.585 AM
% $Revision: 0.1 $  $Date: Tue 12/13/2016 11:58:46.592 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
narginchk(3, 3);

Q = size(P, 1); % number of chirplets

fs = 1; % sampling frequency set to 1 Hz

N0 = n(1);
N = n(end) + 1;

d = (f(end) - f(1))/length(f);
nfreq = round(fs / d);
fr = round(f * nfreq / fs + 1); % frequency range

% =========================================================================
% estimate WVD
% =========================================================================
wig = zeros(length(fr), 2*(N - N0)); 		% bug, should be 2*nfreq?
for k = 1:Q
    chirplet_k = make_chirplets(N, P(k,:));
    tp = real_wigner1(chirplet_k, fs, nfreq);
    
    tp_k = tp(fr, (2*N0 + 1):2*N);
    wig = wig+tp_k;
end %for

end % function chrpltwvd_direct

% [EOF]
