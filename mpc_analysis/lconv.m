function y = lconv(x, h)
% LCONV -- perform a linear convolution with ffts
%
%  Syntax:
%    y = lconv(x, h)
%
%  Inputs:
%    x, h   input vectors
%
%  Outputs:
%    y      the linear convolution of x and h

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; 02-Mar-2005
% $Revision: 0.1 $  $Date: 02-Mar-2005 22:25:16$

narginchk(2, 2);

x = x(:);
N = length(x);
h = h(:);
M = length(h);
P = 2^nextpow2(N+M-1);

y = ifft( fft(x,P) .* fft(h,P));
y = y(1:N+M-1);

if (isreal(x) & isreal(h))
  y = real(y);
end
