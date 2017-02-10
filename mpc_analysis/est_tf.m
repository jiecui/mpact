function [t, f] = est_tf(x, c, d, M, decf)
% EST_TF -- Estimate the location in time and frequency of the chirp.
%
%  Syntax:
%    [t f] = est_tf(x, c, d, M, decf)
%
%  Inputs:
%    x     signal vector
%    c     current estimate of the chirp rate (optional, default is 0)
%    d     current estimate of the duration (optional, default is 5)
%    M     number of points to compute in frequency (optional, default is 64)
%    decf  time decimation factor of the spectrogram (optional, default is 1)
%
%  Outputs:
%    t     estimate of the location in time
%    f     estimate of the location in frequency
%
% Given the current estimates of the chirp rate and duration, estimate
% the location in time and frequency from the maximum value of the
% spectrogram.  M and decf serve to decrease computations but could 
% decrease performance.  Algorithm is O(NM/decf log M)

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; 02-Mar-2005
% $Revision: 0.1 $  $Date: 02-Mar-2005 22:25:16$

narginchk(1, 5);

if (nargin < 2)
  c = 0;
end
if (nargin < 3)
  d = 5;
end
if (nargin < 4)
  M = 64;
end
if (nargin < 5)
  decf = 1;
end

M = 2*floor(M/2);  % want M to be even
x = x(:);
N = length(x);

% compute spectrogram window and prune to save computations
hN = 2*round(32*d/5)-1;  % lessens  truncation of the window
hN = min(hN, 4*M+1);
hN = min(hN, 2*floor(N/2)-1);
h = conj(make_chirplets(hN, [1 (hN+1)/2 0 c d]));

% compute spectrogram and find the location of the maximum
S = spec2(x,1,M,decf,h);
[~, t] = max(max(S,[],1));
[~, f] = max(S(:,t));

% convert from sample number to real units
t = (t-1)*decf + 1;
f = 2*pi*(f-1)/M - pi;
