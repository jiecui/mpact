function d = est_d(x, t, f, c, M, dlow, dhigh)
% EST_D -- Estimate the duration.
%
%  Syntax:
%    d = est_d(x, t, f, c, M, dlow, dhigh)
%
%  Inputs:
%    x      signal vector
%    t      current estimate of the location in time
%    f      current estimate of the location in frequency
%    c      current estimate of the chirp rate
%    M      number of points in the grid (optional, default is 64)
%    dlow   lowest value of the duration (optional, default is 0.25)
%    dhigh  highest value of the duration (optional, default is M/2)
%
%  Outputs
%    d      duration
%
% Given the current estimates of the location in time and frequency and
% chirp rate, estimate the duration using a grid search on the 
% likelihood function.

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; 02-Mar-2005
% $Revision: 0.1 $  $Date: 02-Mar-2005 22:25:16$

x = x(:);
N = length(x);
 
narginchk(4, 7);
if (nargin < 5)
  M = 64;
end
if (nargin < 6)
  dlow = 0.25;
end
if (nargin < 7)
  dhigh = M/2;      % dhigh = M/2
end

% center in time and window, window length = M
% should probably use a non-rectangular window -> gaussian, hamming etc.
MM = M/2;
rt = round(N/2);
if M > N
    xx = zeros(M, 1);
    t0 = MM - rt;
    t1 = t0 + N -1;
    xx(t0:t1) = x;
    x = xx;
elseif M < N
    t0 = rt - MM;
    t1 = t0 + M - 1;
    x = x(t0:t1);
end % if

% create functions corresponding to all possible values of d
y = zeros(M,M);
dd = linspace(dlow, dhigh, M);
for i=1:M
  y(:,i) = make_chirplets(M,[1 M/2 f c dd(i)]);
  y(:,i) = y(:,i)/norm(y(:,i));
end

% find the d that gives the maximum of the likelihood function
A = abs(x'*y);
[~, i] = max(A);
d = dd(i);
