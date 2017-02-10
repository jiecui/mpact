function [freq,wd,amp,phi,tc,fc,cr,d] = instchrppara(P,n,f)
% INSTCHRPPARA -- Get the parameters of the instant prime chirplet
%
% Syntax:
%   [freq,wd,amp,phi,tc,fc,cr,d] = instchrppara(P,n,f)
%
% Inputs:
%   P       matrix of parameters [amp time freq chirp_rate duration; ...]
%   n       vector of time range (assume fs = 1 Hz)
%   f       vector of frequency range ( assume [0 1] range) at each time
%           point
%
% Outputs:
%   freq    frequency position of the max wigner density
%   wd      the max Wigner-Ville distribution density
%   amp     amplitude of the prime chirplet
%   phi     the phase
%   d       duration
%   cr      chirp rate
%   tc      time center
%   fc      frequency center
%
% Assume sampling frequency is 1 Hz.
% ToDo:
%   * to change to arbitrary sampling frequency

% Copyright 2005-2016 Richard J. Cui. Created: Sun 03/13/2005 10:02:34.349 PM
% $Revision: 0.5 $  $Date: Tue 12/13/2016  1:13:50.617 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% parse inputs
narginchk(1,3);

% parameter settings
n = n(:)';          % force to be row vector
N = length(n);
freq = zeros(N,1);
wd = freq; amp = freq; phi = freq;
d = freq; cr = freq; tc = freq; fc = freq;

p = size(P,1);      % number of chirplets

% find the prime parameters
mw = zeros(1,p);
idx = zeros(1,p);
for l = 1:N
    for k = 1:p
        wig = chrpltwvd_explicit(P(k,:),n(l),f);
        [mw(k),idx(k)] = max(wig);
    end %for
    [wd(l),midx] = max(mw);
    freq(l) = f(idx(midx));
    prmP = P(midx,:);       % prime P at this instant n(l)
    amp(l) = abs(prmP(1));
    phi(l) = -angle(prmP(1));
    tc(l) = prmP(2);
    fc(l) = prmP(3);
    cr(l) = prmP(4);
    d(l) = prmP(5);
end %for

end