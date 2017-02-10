function wig = chrpltwvd_explicit(P, n, f)
% CHRPLTWVD_EXPLICIT compute chirplet spectrogram with explicit formulae
%
% Sytax:
%   wig = chrplt_wigner(P, n, f)
%
% Inputs:
%   P       - matrix of parameters [amp time freq chirp_rate duration; ...]
%   n       - vector of time range (assume fs = 1 Hz)
%   f       - vector of frequency range ([0 1] range)
%
% Outputs:
%   wig     - length(n)_by_length(f) WVD
%
% Note:
%   Assume the sampling frequency of chirplet is normalized to 1 Hz.

% Copyright 2005-2016 Richard J. Cui. Created: Sun 03/13/2005 10:02:34.349 PM
% $Revision: 0.5 $  $Date: Tue 12/13/2016  1:13:50.617 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% parse inputs
narginchk(3, 3);

% parameter settings
n = n(:)';      % force to be row vector
f = f(:);       % force to be column vector
% convert to rad in [0 2*pi] range
f = f * 2 * pi;

m = length(f);  % number of freq points
l = length(n);  % number of time points
nchirp = size(P,1);    % number of chirplets

wig = zeros(m, l);
nmat = repmat(n,m,1);  % time matrix: m by l
fmat = repmat(f,1,l);  % freq matrix: m by l
for k = 1:nchirp
    A_k     = P(k,1); % the complex amplitude
    tc_k    = P(k,2); % time center
    fc_k    = P(k,3); % frequency center (rad)
    c_k     = P(k,4); % the chirprate
    d_k     = P(k,5); % duration of the chirplet

    % WVD = w1*w2*w3
    % --------------
    %   w1 = a^2/pi
    %   w2 = exp(-(n-tc)^2/2/d^2)
    %   w3 = exp(-2*d^2*((f-fc)-c*(n-tc))^2)
    % use matrix method. cannot be very large matrix
    a = abs(A_k);
    w1 = a^2/pi;
    w2 = exp(-(nmat-tc_k).^2/2/d_k^2);
    w3 = exp(-2*d_k^2*((fmat-fc_k)-c_k*(nmat-tc_k)).^2);
    w_k = w1*w2.*w3;

    wig = wig + w_k;
end % for

end % funciton

% [EOF]
