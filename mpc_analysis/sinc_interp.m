function out = sinc_interp(x, a)
% SINC_INTERP -- sinc interpolate a signal
%
%  Syntax
%    y = sinc_interp(x, a)
%
%  Inputs
%    x      signal vector
%    a      interpolation factor (optional, default is 2)
%
%  Outputs
%    y     interpolated vector.  If N=length(x) then
%          length(y) = a*N-a+1 and y(1:a:end) = x.
%          No extrapolation is done.
%
% Surprisingly this does not seem to be included with matlab.

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; 02-Mar-2005
% $Revision: 0.1 $  $Date: 02-Mar-2005 22:25:16$

narginchk(1, 2);

if (nargin < 2)
  a = 2;
end

x = x(:);
N = length(x);
M = a*N-a+1;

% y has length: a*N-a+1
y = zeros(M,1);
y(1:a:M) = x;

% h has length: 2*(a*N-a-1)+1
h = sinc([-(N-1-1/a):1/a:(N-1-1/a)]');

% out has length 3*(a*N-a)-1
out = lconv(y, h);

% what we want has length: a*N-a+1
out = out(a*N-a:end-a*N+a+1);
